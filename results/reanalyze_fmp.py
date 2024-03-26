"""
reanalyze_fmp.py

Extract dG data for each edge in an .fmp file at the specified time.

"""
import argparse
import logging
import os
from pprint import pformat
import sys

import numpy as np

from schrodinger.application.desmond.measurement import Measurement
from schrodinger.application.scisol.packages.fep import graph
from schrodinger.application.scisol.packages.fep_gui import sid_report
from schrodinger.utils import log


# Configure logging
try:
    # Schrodinger logging
    logger = log.get_output_logger(__name__)
except NameError:
    # Python builtin fallback
    logger = logging.getLogger(__name__)
logger.format = '%(message)s'
logger.level = logging.WARNING


# Used for sanity check of `cpx - sol` vs. `ddg` interpolation
INTERP_DDG_TOL = 0.001


def get_energy_timeseries_arrays(sr, dg_type='forward'):
    """
    Extract time, complex_dg, solvent_dg, and ddg arrays from the SID report.

    """
    assert all(sr.cpx_timestep_list == sr.sol_timestep_list)
    assert dg_type in ['forward', 'reverse']  # TODO: maybe add reverse/sliding
    t = np.array(sr.cpx_timestep_list)
    dg_type_attr_map = {
        'forward': 'delta_g_forward',
        'reverse': 'delta_g_reverse',
    }
    complex_dg = np.array(getattr(sr, f'cpx_{dg_type_attr_map[dg_type]}'))
    solvent_dg = np.array(getattr(sr, f'sol_{dg_type_attr_map[dg_type]}'))
    ddg = complex_dg - solvent_dg
    return t, complex_dg, solvent_dg, ddg


def recalc_fep(g, end_time, relative=False, dg_type='forward'):
    """
    Relculate dG predictions at `end_time` using SID data stored in the graph.

    Parameters
    ----------
    g: schrodinger.application.scisol.packages.fep.graph.Graph
        The graph that will be copied and have its free energies recalculated.
    end_time: int or float
        The time in nanonseconds where the free energies will calculated. The
        complex_dg, solvent_dg and ddg values will be interpolated from the
        available timepoints in the SID report data.

    Return
    ------
    g_new: schrodinger.application.scisol.packages.fep.graph.Graph
        Graph containing the recalculated energies.
    """
    for e in g.edges():
        logger.info(e.short_id_title)

        # Read SID report data
        sr = sid_report.get_edge_data(e)
        logger.debug("got sid report")
        if sr is not None:
            logger.debug("sid report is not None")

            # Get energy timeseries data arrays
            times, cpx_dgs, sol_dgs, ddgs = \
                get_energy_timeseries_arrays(sr, dg_type=dg_type)
            logger.debug(f"times:\n{pformat(times)}\n")
            logger.debug(f"cpx_dgs:\n{pformat(cpx_dgs)}\n")
            logger.debug(f"sol_dgs:\n{pformat(sol_dgs)}\n")
            logger.debug(f"ddgs:\n{pformat(ddgs)}\n")

            # Interpolate values at new timepoint
            new_cpx_dg = np.interp(end_time, times, cpx_dgs)
            new_sol_dg = np.interp(end_time, times, sol_dgs)
            new_ddg = np.interp(end_time, times, ddgs)

            # Sanity check on interpolation
            assert abs(new_cpx_dg - new_sol_dg - new_ddg) <= INTERP_DDG_TOL

            # Set new dG values
            e.complex_dg.val = new_cpx_dg
            e.solvent_dg.val = new_sol_dg

    # Run cycle closure to generate new bennett_ddg, pred_dg, etc.
    logger.info('calculating cycle closure')
    g.calc_cycle_closure(ignore_exp_dg=relative)

    return g


def main(fmp, time, out_fmp=None, relative=False, dg_type='forward'):
    '''Main workflow'''
    logger.info(f'loading graph {fmp}')
    g = graph.Graph.deserialize(fmp)
    logger.info(f'recalculating energies at {time}ns')
    gg = recalc_fep(g, time, relative, dg_type)
    if out_fmp is None:
        # TODO use os.path.splitext
        base, ext = os.path.splitext(fmp)
        dg_type_suffix = {
            'forward': '',
            'reverse': 'rev',
        }
        out_fmp = f'{base}_reanalyzed_{time}ns{dg_type_suffix[dg_type]}{ext}'
    logger.info(f'writing {out_fmp}')
    gg.write(out_fmp)
    logger.info('done.')


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('fmp',
                        help='input .fmp file'),
    parser.add_argument('-t', '--end-time', required=True, type=int,
                        help='end time for energy recalculation (ns)'),
    parser.add_argument('-o', '--out-fmp', default=None,
                        help='output reanalyzed .fmp filename'),
    parser.add_argument('-d', '--dg_type', default='forward',
                        help='which timeseries value to use '
                             '(currently either "forward" or "reverse")')
    parser.add_argument('--relative', action='store_true',
                        help='set pred_dg values relative to reference node')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('--debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    # Adjust logging
    if args.debug:
        logger.level = logging.DEBUG
    elif args.verbose:
        logger.level = logging.INFO

    main(args.fmp, args.end_time, out_fmp=args.out_fmp, relative=args.relative,
         dg_type=args.dg_type)
