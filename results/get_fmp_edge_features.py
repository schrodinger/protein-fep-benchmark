import argparse
from collections import defaultdict
import logging
import os
import sys

import pandas as pd

from schrodinger import structure
from schrodinger.application.prime.packages.ifd_plus_stages import stage
from schrodinger.application.scisol.fep import graph
from schrodinger.application.scisol.fep.groups.protein_corrections \
    import get_edge_mutated_sites
from schrodinger.infra import mm
from schrodinger.structutils import analyze
from schrodinger.structutils.interactions import hbond
from schrodinger.utils import log

# Local modules
from fmp_features import calculate_graph_features
import mutation


# Configure logging
DEFAULT_LOGGING_LEVEL = logging.WARNING
logger = log.get_output_logger(__name__)
logger.level = DEFAULT_LOGGING_LEVEL


# H-bond params same as for IFD+
MAX_HBOND_DIST = stage.Stage.MAX_HBOND_DIST  # 2.85
MIN_DONOR_ANGLE = stage.Stage.MIN_DONOR_ANGLE  # 90
MIN_ACCEPTOR_ANGLE = stage.Stage.MIN_ACCEPTOR_ANGLE  # 60

# Rescode lists
NEGATIVE_RESIDUES = ['ASP', 'GLU']
POSITIVE_RESIDUES = ['HIP', 'LYS', 'ARG']
DEFAULT_CHARGED_NEIGHBOR_DISTANCES = [5.0, 7.5, 10.0]
DEFAULT_HA_NEIGHBOR_DISTANCE = 5.0
TITRATABLE_RESCODES = "ARG,ASH,ASP,GLH,GLU,HID,HIE,HIP,HIS,LYN,LYS"

# Interpretable secondary structure values
SECONDARY_STRUCTURE_VALUE_TO_STR_MAP = {
    mm.MMCT_SS_NONE: 'none',
    mm.MMCT_SS_LOOP: 'loop',
    mm.MMCT_SS_HELIX: 'helix',
    mm.MMCT_SS_STRAND: 'strand',
    mm.MMCT_SS_TURN: 'turn',
}


# Result dict keys
SB_KEY = 'salt_bridge'
N_SB_KEY = 'n_salt_bridges'
XINT_SB_KEY = 'cross_int_salt_bridge'
N_CHG_NBR_KEY = 'n_charged_neighbors'
N_CPX_PI_HB_KEY = 'n_cpx_pi_hbonds'
N_SOL_PI_HB_KEY = 'n_sol_pi_hbonds'
RES_SS_KEY = 'res_ss'
RES_SS_STR_KEY = 'res_ss_str'
N_CPX_HA_NBR_KEY = "n_cpx_ha_neighbors"
N_SOL_HA_NBR_KEY = "n_sol_ha_neighbors"

PRIME_ENERGY_PROPERTY = 'r_psp_Prime_Energy'


##################################################
# Utility functions
##################################################

def residue_str(res):
    '''Return a standard string representation of the residue object.'''
    return f'{res.chain}:{res.pdbres.strip()}{res.inscode.strip()}'


def title_3to1(t):
    '''
    Return the single-letter version of the node title.
    '''
    pmn = mutation.ProteinMutationNode(t)
    mut1_str = ",".join(pmn.mutations1).replace("-", ":")
    return mut1_str if mut1_str else "WT"


##################################################
# Structural feature functions
##################################################


def get_hbonds(st, res):
    '''Return a list of hbonds.'''
    res_str = f"{res.chain}:{res.pdbres.strip()}{res.resnum}{res.inscode}"
    logger.debug(f'detecting salt bridge for {res_str}')

    res_sidechain_asl = f'({res.getAsl()} and sidechain)'
    res_sidechain_atoms = analyze.evaluate_asl(st, res_sidechain_asl)

    st_atoms = analyze.evaluate_asl(st, 'protein')

    return hbond.get_hydrogen_bonds(st, atoms1=res_sidechain_atoms,
                                      atoms2=st_atoms,
                                      max_dist=MAX_HBOND_DIST,
                                      min_donor_angle=MIN_DONOR_ANGLE,
                                      min_acceptor_angle=MIN_ACCEPTOR_ANGLE,
                                      honor_pbc=False)


def find_salt_bridges(st, res):
    '''
    Return a list of `[donor, acceptor]` H-bond pairs.

    Returned `donor` and `acceptor` are _structure._StructureAtom instances.
    '''
    hbonds = get_hbonds(st, res)

    # res_sidechain_asl = f'({res.getAsl()} and sidechain)'
    # res_sidechain_atoms = analyze.evaluate_asl(st, res_sidechain_asl)

    sb_list = []
    for donor, acceptor in hbonds:
        donor_res = donor.pdbres.strip()
        acceptor_res = acceptor.pdbres.strip()

        donor_str = f'{residue_str(donor)}.{donor.pdbname.strip()}'
        acceptor_str = f'{residue_str(acceptor)}.{acceptor.pdbname.strip()}'
        if donor.index in res.atom:
            res_is_donor = True
            logger.debug(f'- residue is donor: {donor_str}')
        else:
            res_is_donor = False
            logger.debug(f'- residue is acceptor: {acceptor_str}')

        hbond_is_salt_bridge = (donor_res in POSITIVE_RESIDUES) and \
                (acceptor_res in NEGATIVE_RESIDUES)

        if hbond_is_salt_bridge:
            if res_is_donor:
                logger.info('- found salt bridge '
                            f'({donor_str} --- {acceptor_str})')
                sb_list.append([donor, acceptor])
            else:
                logger.info('- found salt bridge '
                            f'({acceptor_str} --- {donor_str})')
                sb_list.append([acceptor, donor])

    return sb_list


def count_heavy_atom_neighbors(st, res, dist=DEFAULT_HA_NEIGHBOR_DISTANCE):
    '''
    Return a count of the number of solute heavy atoms neighbors.

    Count number of non-hydrogen atoms within `dist` Å of sidechain or CA atom,
    excluding the atoms of the residue.
    '''
    res_asl = (f'c. {res.chain} and res.num {res.resnum} and '
               f'''res.ins \"{res.inscode or ' '}\"''')
    neighbor_asl = (
        '(protein and not at.elem H)'
        f' and (within {dist} ({res_asl} and (sidechain or at.ptype CA)))'
        f' and not ({res_asl})'
    )
    neighbor_atoms = analyze.evaluate_asl(st, neighbor_asl)
    n_neighbors = len(neighbor_atoms)
    print(f'→ found {n_neighbors} heavy atom neighbors for '
          f'{res.chain}:{res.pdbres.strip()}{res.resnum}{res.inscode.strip()}')
    return n_neighbors


def get_charged_neighbor_counts(st, res, dist=DEFAULT_CHARGED_NEIGHBOR_DISTANCES):
    '''
    Return the number of charged/titratable neighbors for the residue for each
    distance cutoff in the `dist` list.
    '''
    res_str = f'{res.chain}:{res.pdbres.strip()}{res.resnum}{res.inscode.strip()}'
    logger.debug(f'getting # of charged neighbors for {res_str}')
    neighbor_counts = {}
    res_asl = f"c. {res.chain} and res.num {res.resnum} and res.ins \"{res.inscode or ' '}\""
    for d in dist:
        nearby_sidechains = f'(sidechain) and (within {d} ({res_asl} and sidechain))'
        titr = f'res.ptype {TITRATABLE_RESCODES}'
        asl = f"(fillres (({nearby_sidechains}) and ({titr}))) and at.ptype CA"
        atoms = analyze.evaluate_asl(st, asl)
        n_neighbors = len(atoms)
        neighbor_counts[d] = n_neighbors
        logger.debug(f'found {n_neighbors} charged neighbors within {d} Å')
    return neighbor_counts





def find_pi_hbonds(st, res, max_dist=3):
    '''
    Return the number of potential pi Hbonds for the residue.
    '''
    res_str = f'{res.chain}:{res.pdbres.strip()}{res.resnum}{res.inscode.strip()}'
    logger.debug(f'getting # of pi Hbond partners for {res_str}')

    sc_asl = f'({res.getAsl()}) and sidechain' #and not (at.elem H) and smarts. [R])
    ring_asl = 'smarts. [R] and not at.elem H'
    polar_h_asl = 'at.elem H and not (/C0-H0/)'
    sc_ring_asl = f'({sc_asl}) and ({ring_asl})'
    sc_polar_h_asl = f'({sc_asl}) and ({polar_h_asl})'

    # Potential pi acceptors are ring atoms, not part of the same sidechain,
    # within max_dist cutoff of sidechain polar H.
    nearby_ring_asl = (f'({ring_asl}) and not ({sc_asl}) '
                       f'and (within {max_dist} ({sc_polar_h_asl}))')

    # Potential donors are polar H, not part of the same sidechain,
    # and within max_dist_cutoff of the sidechain ring atoms.
    nearby_polar_h_asl = (f'({polar_h_asl}) and not ({sc_asl}) '
                          f'and (within {max_dist} ({sc_ring_asl}))')

    # Get the atom id lists
    sc_ring_atom_ids = analyze.evaluate_asl(st, sc_ring_asl)
    sc_polar_h_atom_ids = analyze.evaluate_asl(st, sc_polar_h_asl)
    nearby_ring_atom_ids = analyze.evaluate_asl(st, nearby_ring_asl)
    nearby_polar_h_atom_ids = analyze.evaluate_asl(st, nearby_polar_h_asl)

    n_pi_hbonds = 0

    # Polar Hbond donors to sidechain ring
    if len(sc_ring_atom_ids):
        print(f'found sidechain ring for {res_str}')
        print(f'- nearby polar H atom indices: {nearby_polar_h_atom_ids}')
        nearby_polar_h_residues = list(set([
            st.atom[i].getResidue()
            for i in nearby_polar_h_atom_ids
        ]))
        n_pi_hbonds += len(nearby_polar_h_residues)

    # Aromatic pi acceptors of Hbond from sidechain polar H
    if len(sc_polar_h_atom_ids):
        print(f'found sidechain polar H for {res_str}')
        print(f'- nearby ring atom indices: {nearby_ring_atom_ids}')
        nearby_ring_residues = list(set([
            st.atom[i].getResidue()
            for i in nearby_ring_atom_ids
        ]))
        n_pi_hbonds += len(nearby_ring_residues)

    if n_pi_hbonds:
        print(f'- found {n_pi_hbonds} pi Hbonds!')

    return n_pi_hbonds


def process_edge_position(e, pos):
    '''
    Perform analyses for the mutated position along the edge.
    '''
    position_data = []
    for n in e:
        logger.info(f'Node: {n.struc.title}')

        # sol/cpx structures and residue objects
        sol_st = n.struc
        cpx_st = n.graph.receptor_struc.merge(n.struc)
        sol_res = sol_st.findResidue(pos)
        cpx_res = cpx_st.findResidue(pos)

        # detect salt bridges and other Hbond types
        cpx_sb_list = find_salt_bridges(cpx_st, cpx_res)
        sol_sb_list = find_salt_bridges(sol_st, sol_res)

        # cpx_arom_hb_list = find_aromatic_hbonds(cpx_st, cpx_res)
        # sol_arom_hb_list = find_aromatic_hbonds(sol_st, sol_res)

        res_ss = sol_res.secondary_structure
        res_ss_str = SECONDARY_STRUCTURE_VALUE_TO_STR_MAP[res_ss]

        n_cpx_ha_neighbors = count_heavy_atom_neighbors(cpx_st, cpx_res)
        n_sol_ha_neighbors = count_heavy_atom_neighbors(sol_st, sol_res)

        node_result = {
            SB_KEY: len(cpx_sb_list) > 0,
            N_SB_KEY: len(cpx_sb_list),
            XINT_SB_KEY: len(cpx_sb_list) > len(sol_sb_list),
            N_CHG_NBR_KEY: get_charged_neighbor_counts(cpx_st, cpx_res),
            N_CPX_PI_HB_KEY: find_pi_hbonds(cpx_st, cpx_res),
            N_SOL_PI_HB_KEY: find_pi_hbonds(sol_st, sol_res),
            RES_SS_KEY: res_ss,
            RES_SS_STR_KEY: res_ss_str,
            N_CPX_HA_NBR_KEY: n_cpx_ha_neighbors,
            N_SOL_HA_NBR_KEY: n_sol_ha_neighbors,
        }

        position_data.append(node_result)
    return position_data


def mutation_breaks_interaction(position_result, interaction_key):
    '''
    Determine whether the mutation breaks the specified action.

    Input format:
    position_result = List[Dict] (i.e. [<n0_result_dict>, <n1_result_dict>])
    interaction_key = str
    '''
    d0, d1 = position_result
    return d0[interaction_key] and not d1[interaction_key]


def mutation_forms_interaction(position_result, interaction_key):
    '''
    Determine whether the mutation forms the specified interaction.

    Same input formats as `mutation_breaks_interaction`
    '''
    d0, d1 = position_result
    return d1[interaction_key] and not d0[interaction_key]


def get_max_charged_neighbors_by_cutoff(edge_result, key='n_charged_neighbors'):
    '''
    Return the max # of charged neighbors
    '''
    cutoff_n_neighbor_lists = defaultdict(list)
    for position_result in edge_result:
        for node_dict in position_result:
            for d, n in node_dict[key].items():
                cutoff_n_neighbor_lists[d].append(n)
    return {d: max(n) for d, n in cutoff_n_neighbor_lists.items()}


##################################################
# Edge/node feature functions
##################################################

def get_prime_energy(n):
    '''Return the "Prime Energy" property for the Node.'''
    return n.struc.property.get(PRIME_ENERGY_PROPERTY)


def residue_type(res):
    '''return the residue type for the residue'''
    three_letter_resname = res.pdbres.strip()
    rescode = structure.RESIDUE_MAP_3_TO_1_LETTER[three_letter_resname]
    if rescode in 'AVILM':
        return 'nonpolar'
    elif rescode in 'CSTNQ':
        return 'polar'
    elif rescode in 'DE':
        return 'negative'
    elif rescode in 'HKR':
        return 'positive'
    elif rescode in 'FYW':
        return 'aromatic'
    elif rescode == 'G':
        return 'glycine'
    elif rescode == 'P':
        return 'proline'
    else:
        raise ValueError('unrecognized rescode')


##################################################
# Workflow
##################################################

def consolidate_edge_position_results(edge_position_results_list):
    '''
    Analyze and consolidate position results for the edge into a single dict.
    suitable for writing to CSV.

    '''
    breaks_salt_bridge = any(
        mutation_breaks_interaction(position_result, SB_KEY)
        for position_result in edge_position_results_list
    )
    forms_salt_bridge = any(
        mutation_forms_interaction(position_result, SB_KEY)
        for position_result in edge_position_results_list
    )
    breaks_cross_int_salt_bridge = any(
        mutation_breaks_interaction(position_result, XINT_SB_KEY)
        for position_result in edge_position_results_list
    )
    forms_cross_int_salt_bridge = any(
        mutation_forms_interaction(position_result, XINT_SB_KEY)
        for position_result in edge_position_results_list
    )
    res_ss = ','.join(list(set([
        str(p[RES_SS_KEY])
        for position_result in edge_position_results_list
        for p in position_result
    ])))
    res_ss_str = ','.join(list(set([
        str(p[RES_SS_STR_KEY])
        for position_result in edge_position_results_list
        for p in position_result
    ])))

    start_n_cpx_pi_hbonds, end_n_cpx_pi_hbonds = [
        p[N_CPX_PI_HB_KEY]
        for position_result in edge_position_results_list
        for p in position_result
    ]
    start_n_sol_pi_hbonds, end_n_sol_pi_hbonds = [
        p[N_SOL_PI_HB_KEY]
        for position_result in edge_position_results_list
        for p in position_result
    ]

    start_n_cpx_ha_neighbors, end_n_cpx_ha_neighbors = [
        p[N_CPX_HA_NBR_KEY]
        for position_result in edge_position_results_list
        for p in position_result
    ]
    start_n_sol_ha_neighbors, end_n_sol_ha_neighbors = [
        p[N_SOL_HA_NBR_KEY]
        for position_result in edge_position_results_list
        for p in position_result
    ]

    consolidated_results = {
        'breaks_salt_bridge': breaks_salt_bridge,
        'forms_salt_bridge': forms_salt_bridge,
        'breaks_cross_int_salt_bridge': breaks_cross_int_salt_bridge,
        'forms_cross_int_salt_bridge': forms_cross_int_salt_bridge,
        RES_SS_KEY: res_ss,
        RES_SS_STR_KEY: res_ss_str,
        'start_n_cpx_pi_hbonds': start_n_cpx_pi_hbonds,
        'end_n_cpx_pi_hbonds': end_n_cpx_pi_hbonds,
        'start_n_sol_pi_hbonds': start_n_sol_pi_hbonds,
        'end_n_sol_pi_hbonds': end_n_sol_pi_hbonds,
        f'start_{N_CPX_HA_NBR_KEY}': start_n_cpx_ha_neighbors,
        f'end_{N_CPX_HA_NBR_KEY}': end_n_cpx_ha_neighbors,
        f'start_{N_SOL_HA_NBR_KEY}': start_n_sol_ha_neighbors,
        f'end_{N_SOL_HA_NBR_KEY}': end_n_sol_ha_neighbors,
    }

    # Get max number of charged neighbors at each distance cutoff
    charged_neighbors = get_max_charged_neighbors_by_cutoff(
        edge_position_results_list
    )
    for d, n in charged_neighbors.items():
        consolidated_results[f'charged_neighbors_{d}'] = n

    return consolidated_results


def process_edge(e):
    '''
    Perform all analyses for the edge.
    '''
    logger.info('\n' + '=' * 50)
    logger.info(f'Edge {e.short_id_title}')
    logger.info('=' * 50)
    positions = [site.as_str() for site in get_edge_mutated_sites(e)]
    print(e.short_id_title, positions)
    logger.debug(f'positions = {positions}')

    n0, n1 = e

    # Get position results
    edge_position_results = []
    for p in positions:
        position_result = process_edge_position(e, p)
        edge_position_results.append(position_result)

    # Prime energy
    delta_prime_energy = get_prime_energy(n1) - get_prime_energy(n0)

    # Assemble edge and position data
    edge_result = {
        'n0': n0.struc.title,
        'n1': n1.struc.title,
        'n0_aa1': title_3to1(n0.struc.title),
        'n1_aa1': title_3to1(n1.struc.title),
        'delta_prime_energy': delta_prime_energy,

        **consolidate_edge_position_results(edge_position_results)
    }

    return edge_result


def main(fmp, write=True, force=False):
    '''Main workflow'''
    # Configure output file
    fmp_basename = os.path.basename(fmp)
    system = '_'.join(fmp_basename.split('_')[:2])  # brittle but functional
    csv_out = f'features/features_{system}.csv'

    # Don't overwrite existing files
    if write and os.path.isfile(csv_out) and not force:
        logger.error(f"Output file {csv_out} exists!  Exiting...")
        sys.exit()

    g = graph.Graph.deserialize(fmp)

    # Limit edge data to direct mutations from WT
    edge_results = [
        process_edge(e)
        for e in g.edges_iter()
        if e[0].struc.title.startswith("WT")
    ]

    df1 = pd.DataFrame.from_dict(edge_results)

    # Create a `system` column at the front
    df1.insert(0, 'system', system)

    # Gather extended fmp features
    ext_csv_out = f'features/ext_features_{system}.csv'

    # Use existing files when possible - extended features are more expensive
    # to calculate.
    if os.path.isfile(ext_csv_out) and not force:
        logger.warning(f"Using existing output file {ext_csv_out}")
        df2 = pd.read_csv(ext_csv_out)
    else:
        logger.warning(f"Calculating extended fmp features for {system}...")
        extended_features = calculate_graph_features(g, logger)
        df2 = pd.DataFrame(extended_features)
        df2.insert(0, 'system', system)
        if write:
            logger.warning(f"Writing extended fmp features to {ext_csv_out}")
            df2.to_csv(ext_csv_out, index=False)
        else:
            logger.warning(
                f"DRY RUN - Would write extended fmp features to {ext_csv_out}"
            )

    df = df2.merge(df1, on=['system', 'n0', 'n1' ], how='left')
    logger.info(df)
    if write:
        logger.warning(f"Writing full fmp features to {csv_out}")
        df.to_csv(csv_out, index=False)
    else:
        logger.warning(f"DRY RUN - Would write full fmp features to {csv_out}")
    logger.warning(f"Done")


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('fmp',
                        help='input .fmp file to analyze')
    parser.add_argument('-n', '--dry-run', dest="write", action="store_false",
                        help="don't write output files (overrides -f/--force)")
    parser.add_argument('-f', '--force', action="store_true",
                        help="overwrite existing output files")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('--debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    main(args.fmp, write=args.write, force=args.force)
