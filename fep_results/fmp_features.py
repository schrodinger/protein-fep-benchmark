'''
fmp_features.py

Library for extracting feature information from an .fmp file.
'''
import json
from pprint import pformat
import statistics

import numpy as np
import scipy

from schrodinger.application.desmond import new_fep_edge_data as fed
from schrodinger.application.bioluminate.protein import PropertyCalculator
from schrodinger.application.desmond.fep_edge_report_maker import parse_res_tag
from schrodinger.application.scisol.packages.fep_gui import sid_report
from schrodinger.structutils import analyze

# Custom local module (depends on SCHRODINGER_PYTHONPATH)
from . import mutation


# Theoretical Maximum SASA values from Table 1 of:
# Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J. & Wilke, C. O.
# "Maximum Allowed Solvent Accessibilites of Residues in Proteins."
# Plos One 8, e80635 (2013).
# https://doi.org/10.1371/journal.pone.0080635
RESIDUE_MAX_SASA_THEORETICAL_TIEN_2013 = {
    'ALA': 129.0,
    'ARG': 274.0,
    'ASN': 195.0,
    'ASP': 193.0,
    'ASH': 193.0,
    'CYS': 167.0,
    'GLU': 223.0,
    'GLH': 223.0,
    'GLN': 225.0,
    'GLY': 104.0,
    'HID': 224.0,
    'HIE': 224.0,
    'HIP': 224.0,
    'HIS': 224.0,
    'ILE': 197.0,
    'LEU': 201.0,
    'LYS': 236.0,
    'MET': 224.0,
    'PHE': 240.0,
    'PRO': 159.0,
    'SER': 155.0,
    'THR': 172.0,
    'TRP': 285.0,
    'TYR': 263.0,
    'VAL': 174.0,
}
# Miller values are used by SKEMPI database
RESIDUE_MAX_SASA_MILLER_1987 = {
    'ALA': 113.0,
    'ARG': 241.0,
    'ASN': 158.0,
    'ASP': 151.0,
    'ASH': 151.0,
    'CYS': 140.0,
    'GLU': 183.0,
    'GLH': 183.0,
    'GLN': 189.0,
    'GLY': 85.0,
    'HID': 194.0,
    'HIE': 194.0,
    'HIP': 194.0,
    'HIS': 194.0,
    'ILE': 182.0,
    'LEU': 180.0,
    'LYN': 211.0,
    'LYS': 211.0,
    'MET': 204.0,
    'PHE': 218.0,
    'PRO': 143.0,
    'SER': 122.0,
    'THR': 146.0,
    'TRP': 259.0,
    'TYR': 229.0,
    'VAL': 160.0,
}
RESIDUE_MAX_SASA = RESIDUE_MAX_SASA_MILLER_1987

MIN_SLOPE_TIME_NS = 5
DEFAULT_SLOPE_FRAC = 0.25

DG_PRECISION = 2
LQMS_PRECISION = 4
SASA_PRECISION = 4
INT_FREQ_PRECISION = 4

SID_INTERACTION_TYPES = [
    'hbond', 'l_hbond', 'p_metal', 'l_metal', 'hphob',
    'pipi', 'picat', 'ionic', 'wat_br', 'lig_wat',
]


def calc_sasa_fsasa(res, prop_calc=None):
    '''
    Calculate residue SASA and fractional SASA.
    '''
    st = res.structure
    if prop_calc is None:
        prop_calc = PropertyCalculator(st, 'tmp')
    assert st == prop_calc.struct
    aa = res.pdbres.strip()
    sasa = prop_calc.getResidueSASA(res, sidechain=False)
    sc_sasa = prop_calc.getResidueSASA(res, sidechain=True)
    try:
        fsasa = sasa / RESIDUE_MAX_SASA[aa]
        sc_fsasa = sc_sasa / RESIDUE_MAX_SASA[aa]
    except KeyError:
        # noncanonical amino acid
        fsasa = None
        sc_fsasa = None
    return sasa, fsasa, sc_sasa, sc_fsasa


def determine_structure_location(sol_fsasa, cpx_fsasa, cutoff=1e-4):
    '''
    Determine the structure location for a residue based on fSASA values.
    '''
    if None in [cpx_fsasa, sol_fsasa]:
        # noncanonical AAs have fsasa = None
        return None

    delta_fsasa = cpx_fsasa - sol_fsasa

    if abs(delta_fsasa) < cutoff:
        # no change in fsasa, not an interface position
        if cpx_fsasa < 0.25:
            # buried
            loc = "INT"
        else:
            # solvent-exposed
            loc = "SUR"
    else:
        # interface position
        if sol_fsasa < 0.25:
            # already (mostly) buried in solvent
            loc = "SUP"
        elif cpx_fsasa < 0.25:
            # only buried in complex
            loc = "COR"
        else:
            # still mostly solvent exposed in complex
            loc = "RIM"
    return loc


# def get_residue_size(rescode):
#     # For some reason the res codes in structure._res_sizes are 4 characters
#     if len(rescode) == 3:
#         rescode = rescode + ' '
#     try:
#         return structure._res_sizes[rescode]
#     except:
#         msg = f'Only standard residue codes are supported, skipping {rescode}'
#         logger.warning(msg)
#         return None


def get_energy_timeseries(sr, dg_type='forward'):
    assert all(sr.cpx_timestep_list == sr.sol_timestep_list)
    assert dg_type in ['forward']  # TODO: maybe add reverse/sliding
    t = np.array(sr.cpx_timestep_list)
    dg_type_attr_map = {
        'forward': 'delta_g_forward',
    }
    complex_dg = np.array(getattr(sr, f'cpx_{dg_type_attr_map[dg_type]}'))
    solvent_dg = np.array(getattr(sr, f'sol_{dg_type_attr_map[dg_type]}'))
    ddg = complex_dg - solvent_dg
    # return pd.DataFrame(
    #     list(zip(t, complex_dg, solvent_dg, ddg)),
    #     columns=['time_ns', 'complex_dg', 'solvent_dg', 'ddg']
    # )
    return t, complex_dg, solvent_dg, ddg


def get_pipi_interaction_dict_from_stats(stats, precision=4):
    '''
    `stats` for pi-pi interactions looks like this:
    (
        0.068,                           # frequency / persistence
        'B:TYR_14',                      # protein residue label
        [4.0, -22.1, -8.7],              # Calpha xyz coords
        [3698, 3699, 3701, 3702, 3700],  # ring atoms on ligand (mutated res)
        0.0                              # face-to-face fraction
    )
    '''
    freq, res_label, _, _, f2f_frac = stats
    return {
        'frequency': float(freq),
        'interacting_residue': dict(zip(
            ['chain', 'resname', 'full_resnum', 'atom'],
            parse_res_tag(res_label),
        )),
        'face_to_face_fraction': float(f2f_frac),
    }


def get_picat_interaction_dict_from_stats(stats, precision=4):
    '''
    `stats` for pi-cation interactions looks like one of these:

    Aromatic group on ligand, cation on protein:
     (
        0.16,                      # frequency
        'B:ARG_77:CZ',             # interacting residue label
        [-5.9, -11.1, -5.9],       # interacting residue Calpha coords
        [698, 699, 701, 702, 700]  # ligand pi system atoms
    )

    Cation on ligand, aromatic group on protein:
    (
        1.0,                               # frequency
        'B:TYR_99',                        # interacting residue label
        [13.404797, 3.522047, 15.211519],  # interacting residue Calpha coords
        2945                               # ligand cation atom id
    )
    '''
    freq, res_label, _, aid = stats
    return {
        'frequency': float(freq),
        'interacting_residue': dict(zip(
            ['chain', 'resname', 'full_resnum', 'atom'],
            parse_res_tag(res_label),
        )),
    }


def lq_median_dg_slopes_variances(times, dgs):
    '''
    Calculate last quarter median slope and variance at each timepoint.

    '''
    timestep = times[1] - times[0]
    min_steps = int(np.ceil(MIN_SLOPE_TIME_NS / timestep))
    # at least min_steps, but not longer than the input
    # steps = []
    slopes = []
    variances = []
    for i in range(len(times)):
        if i < 2:
            # n_steps = 0
            slope = np.nan
            var = np.nan
        else:
            cut_times = times[0:i]
            cut_dgs = dgs[0:i]
            n_steps = int(np.ceil(DEFAULT_SLOPE_FRAC * i))
            n_steps = max(n_steps, min_steps)
            n_steps = min(n_steps, i)
            x = cut_times[-n_steps:]
            y = cut_dgs[-n_steps:]
            # Theil-Sen non parametric regression returns a tuple:
            # (median slope, median intercept, slope lower bound, slope upper bound)
            slope = scipy.stats.theilslopes(y, x)[0]
            var = statistics.variance(y)
        # steps.append(n_steps)
        slopes.append(slope)
        variances.append(var)
    return np.array(slopes), np.array(variances)


class NodeFeatures():
    '''
    Class to store feature info from a `graph.Node` instance.
    '''
    def __init__(self, n):
        self._node = n
        self._strucs = {
            'solvent': n.struc,
            'complex': n.graph.receptor_struc.merge(n.struc),
        }
        self._parched_strucs = {}
        # self._position_properties = {}
        self._prop_calc = {}

    def leg_struc(self, leg, parch=False):
        '''
        Return the `structure.Structure` instance for the given leg.

        Optionally remove waters with `parch=True`.
        '''
        if parch:
            try:
                return self._parched_strucs[leg]
            except:
                wet_st = self._strucs[leg]
                water_indices = analyze.evaluate_asl(wet_st, 'water')
                dry_st = wet_st.copy()
                dry_st.deleteAtoms(water_indices)
                assert wet_st.atom_total == dry_st.atom_total + len(water_indices)
                self._parched_strucs[leg] = dry_st
                return dry_st
        return self._strucs[leg]

    def leg_prop_calc(self, leg, parch=False):
        '''
        Return the 'PropertyCalculator' instance for the given leg.
        '''
        try:
            return self._prop_calc[leg]
        except KeyError:
            pc_title = f'{self._node.struc.title}_{leg}'
            pc = PropertyCalculator(self.leg_struc(leg, parch=parch),
                                    pc_title)
            self._prop_calc[leg] = pc
            return pc


def calculate_leg_features_at_position(node_features, leg, pos):
    '''
    Calculate leg-specific structural features for the given position.
    '''
    # Use parched structure for SASA values
    res = node_features.leg_struc(leg, parch=True)
    pc = node_features.leg_prop_calc(leg, parch=True)
    sasa, fsasa, sc_sasa, sc_fsasa = calc_sasa_fsasa(res, pc)
    return {
        f'{leg}_sasa': sasa,
        f'{leg}_fsasa': fsasa,
        f'{leg}_sc_sasa': sc_sasa,
        f'{leg}_sc_fsasa': sc_fsasa,
    }


def calculate_edge_features(e, node_features_map=None, parch=True, logger=None):
    '''
    Calculate sequence- and structure-based features for a single node.

    inputs:
    n = graph.Node instance
    prop_calc (optional) = PropertyCalculator instance

    TODO: remove logging
    '''
    n0, n1 = e

    # Get mutation info
    pme = mutation.ProteinMutationEdge(n0.struc.title, n1.struc.title, e)
    if logger: logger.info(pme)  # (perturbation: {pme.mutation})')

    # Only consider single mutation edges for now
    if not len(pme.mutation) == 1:
        if logger: logger.info('multiple mutations found...skipping')
        return None

    mut = pme.mutation
    pos = mut.position

    node_data = {}
    for n in [n0, n1]:
        prefix = 'start' if n == n0 else 'end'

        # Use existing NodeFeatures if already calculated
        try:
            nf = node_features_map[n]
        except (KeyError, TypeError):
            nf = NodeFeatures(n)

        # for leg in ['complex', 'solvent']:
        #     leg_features = calculate_leg_features(nf, leg)
        #     node_features_dict = {**node_features_dict, **leg_features}

        cpx_res = nf.leg_struc('complex', parch=parch).findResidue(pos)
        cpx_pc = nf.leg_prop_calc('complex', parch=parch)
        cpx_sasa, cpx_fsasa, cpx_sc_sasa, cpx_sc_fsasa = \
            calc_sasa_fsasa(cpx_res, cpx_pc)

        sol_res = nf.leg_struc('solvent', parch=parch).findResidue(pos)
        sol_pc = nf.leg_prop_calc('solvent', parch=parch)
        sol_sasa, sol_fsasa, sol_sc_sasa, sol_sc_fsasa = \
            calc_sasa_fsasa(sol_res, sol_pc)

        delta_sasa = cpx_sasa - sol_sasa
        delta_sc_sasa = cpx_sc_sasa - sol_sc_sasa
        try:
            delta_fsasa = cpx_fsasa - sol_fsasa
            delta_sc_fsasa = cpx_sc_fsasa - sol_sc_fsasa
        except TypeError:
            # noncanonical AAs have fsasa = None (unknown max value)
            delta_fsasa = None
            delta_sc_fsasa = None

        # fSASA-based location (INT, SUR, SUP, RIM, COR)
        loc = determine_structure_location(sol_fsasa, cpx_fsasa)

        # residue size
        size = len([a for a in sol_res.atom if a.element != 'H'])
        formal_charge = sum([a.formal_charge for a in sol_res.atom])

        if logger: logger.debug(f'- {n.struc.title}\n'
                    f'    sasa:   sol = {sol_sasa}\n'
                    f'            cpx = {cpx_sasa}\n'
                    f'          delta = {delta_sasa}\n'
                    f'    fsasa:  sol = {sol_fsasa}\n'
                    f'            cpx = {cpx_fsasa}\n'
                    f'          delta = {delta_fsasa}\n'
                    f'    loc = {loc}')

        node_data = {
            **node_data,
            f'{prefix}_cpx_sasa': round(cpx_sasa, SASA_PRECISION),
            f'{prefix}_cpx_fsasa': None if cpx_fsasa is None else round(cpx_fsasa, SASA_PRECISION),
            f'{prefix}_sol_sasa': round(sol_sasa, SASA_PRECISION),
            f'{prefix}_sol_fsasa': None if sol_fsasa is None else round(sol_fsasa, SASA_PRECISION),
            f'{prefix}_delta_sasa': round(delta_sasa, SASA_PRECISION),
            f'{prefix}_delta_fsasa': None if delta_fsasa is None else round(delta_fsasa, SASA_PRECISION),
            f'{prefix}_cpx_sc_sasa': round(cpx_sc_sasa, SASA_PRECISION),
            f'{prefix}_cpx_sc_fsasa': None if cpx_sc_fsasa is None else round(cpx_sc_fsasa, SASA_PRECISION),
            f'{prefix}_sol_sc_sasa': round(sol_sc_sasa, SASA_PRECISION),
            f'{prefix}_sol_sc_fsasa': None if sol_sc_fsasa is None else round(sol_sc_fsasa, SASA_PRECISION),
            f'{prefix}_delta_sc_sasa': round(delta_sc_sasa, SASA_PRECISION),
            f'{prefix}_delta_sc_fsasa': None if delta_sc_fsasa is None else round(delta_sc_fsasa, SASA_PRECISION),
            f'{prefix}_loc': loc,
            f'{prefix}_size': int(size),
            f'{prefix}_formal_charge': formal_charge,
        }

    # Edge properties

    # Defaults
    missing_array_default = np.array([])
    times      = missing_array_default
    cpx_dgs    = missing_array_default
    sol_dgs    = missing_array_default
    ddgs       = missing_array_default
    cpx_lqms   = missing_array_default
    cpx_lqvars = missing_array_default
    sol_lqms   = missing_array_default
    sol_lqvars = missing_array_default
    ddg_lqms   = missing_array_default
    ddg_lqvars = missing_array_default
    max_wt_pipi_freq = None
    sum_wt_pipi_freq = None
    max_mut_pipi_freq = None
    sum_mut_pipi_freq = None
    max_wt_picat_freq = None
    sum_wt_picat_freq = None
    max_mut_picat_freq = None
    sum_mut_picat_freq = None

    # Get energy timeseries
    sr = sid_report.get_edge_data(e)
    if sr is not None:
        times, cpx_dgs, sol_dgs, ddgs = get_energy_timeseries(sr)

        # Basic max sim time and timestep data
        assert sr.cpx_sim_time == sr.sol_sim_time
        assert sr.cpx_timestep_interval == sr.sol_timestep_interval
        sim_time = sr.cpx_sim_time
        if logger: logger.debug(f'sim time {sim_time}')

        cpx_lqms, cpx_lqvars = lq_median_dg_slopes_variances(times, cpx_dgs)
        sol_lqms, sol_lqvars = lq_median_dg_slopes_variances(times, sol_dgs)
        ddg_lqms, ddg_lqvars = lq_median_dg_slopes_variances(times, ddgs)

        try:
            exp_ddg = pme.graph_edge.exp_ddg.val
        except AttributeError:
            exp_ddg = None

        # Extract interaction stats
        wt_stats_list, mut_stats_list = sr.cpx_sid_lp_results

        # Pi-pi interaction stats
        wt_pipi_dicts = [
            get_pipi_interaction_dict_from_stats(stats)
            for stats in wt_stats_list['pipi']
        ]
        wt_pipi_freqs = [d['frequency'] for d in wt_pipi_dicts]
        # Add a zero term to avoid "empty sequence" ValueError
        max_wt_pipi_freq = max(wt_pipi_freqs + [0])
        sum_wt_pipi_freq = sum(wt_pipi_freqs)

        mut_pipi_dicts = [
            get_pipi_interaction_dict_from_stats(stats)
            for stats in mut_stats_list['pipi']
        ]
        mut_pipi_freqs = [d['frequency'] for d in mut_pipi_dicts]
        max_mut_pipi_freq = max(mut_pipi_freqs + [0])
        sum_mut_pipi_freq = sum(mut_pipi_freqs)

        if logger:
            logger.debug('wt pipi dicts')
            logger.debug(pformat(wt_pipi_dicts, sort_dicts=False))
            logger.debug('mut pipi dicts')
            logger.debug(pformat(mut_pipi_dicts, sort_dicts=False))

        # Pi-cation interaction stats
        wt_picat_dicts = [
            get_picat_interaction_dict_from_stats(stats)
            for stats in wt_stats_list['picat']
        ]
        wt_picat_freqs = [d['frequency'] for d in wt_picat_dicts]
        # Add a zero term to avoid "empty sequence" ValueError
        max_wt_picat_freq = max(wt_picat_freqs + [0])
        sum_wt_picat_freq = sum(wt_picat_freqs)

        mut_picat_dicts = [
            get_picat_interaction_dict_from_stats(stats)
            for stats in mut_stats_list['picat']
        ]
        mut_picat_freqs = [d['frequency'] for d in mut_picat_dicts]
        max_mut_picat_freq = max(mut_picat_freqs + [0])
        sum_mut_picat_freq = sum(mut_picat_freqs)

        if logger:
            logger.debug('wt picat dicts')
            logger.debug(pformat(wt_picat_dicts, sort_dicts=False))
            logger.debug('mut picat dicts')
            logger.debug(pformat(mut_picat_dicts, sort_dicts=False))

        # Ligand (i.e. mutating residue) SASA
        sasa_attrs = [
            'ligand1_cpx_sid_sasa', 'ligand1_sol_sid_sasa',
            'ligand2_cpx_sid_sasa', 'ligand2_sol_sid_sasa',
        ]

    # Create edge data dict
    d = {
        'edge_id': '_'.join(e.short_id),
        'n0': pme.node_titles[0],
        'n1': pme.node_titles[1],
        'mutation': str(mut),
        'exp_ddg': exp_ddg,
        'n_mut': len(mut),
        'start_aa': mut.start,
        'end_aa': mut.end,
        'start_aa1': mut.start1,
        'end_aa1': mut.end1,
        'is_charged': mut.is_charged,
        'is_titration': pme.is_titration,
        'max_wt_pipi_freq': round(max_wt_pipi_freq, INT_FREQ_PRECISION),
        'sum_wt_pipi_freq': round(sum_wt_pipi_freq, INT_FREQ_PRECISION),
        'max_mut_pipi_freq': round(max_mut_pipi_freq, INT_FREQ_PRECISION),
        'sum_mut_pipi_freq': round(sum_mut_pipi_freq, INT_FREQ_PRECISION),
        'max_wt_picat_freq': round(max_wt_picat_freq, INT_FREQ_PRECISION),
        'sum_wt_picat_freq': round(sum_wt_picat_freq, INT_FREQ_PRECISION),
        'max_mut_picat_freq': round(max_mut_picat_freq, INT_FREQ_PRECISION),
        'sum_mut_picat_freq': round(sum_mut_picat_freq, INT_FREQ_PRECISION),
        **node_data,
    }

    # TODO: make this optional
    include_arrays = True
    if include_arrays:
        d = {
            **d,
            'time_arr': json.dumps(times.tolist()),
            'complex_dg_arr': json.dumps(cpx_dgs.tolist()),
            'solvent_dg_arr': json.dumps(sol_dgs.tolist()),
            'ddg_arr': json.dumps(ddgs.tolist()),
            'complex_lqms_arr': json.dumps(cpx_lqms.tolist()),
            'complex_lqvar_arr': json.dumps(cpx_lqvars.tolist()),
            'solvent_lqms_arr': json.dumps(sol_lqms.tolist()),
            'solvent_lqvar_arr': json.dumps(sol_lqvars.tolist()),
            'ddg_lqms_arr': json.dumps(ddg_lqms.tolist()),
            'ddg_lqvar_arr': json.dumps(ddg_lqvars.tolist()),
        }

    return d


def calculate_graph_features(g, logger):
    '''
    Calculate sequence- and structure-based features for the graph.
    '''
    # Create each node struc/prop_calc only once
    # TODO refactor into a class?
    node_features = {
        n: NodeFeatures(n)
        for n in g.nodes()
    }

    g_data = []
    for e in g.edges_iter():
        edge_features = calculate_edge_features(e, node_features, logger=logger)
        if edge_features is not None:
            g_data.append(edge_features)
    return g_data



