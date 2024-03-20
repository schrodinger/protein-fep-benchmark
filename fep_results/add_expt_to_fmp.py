import argparse
import logging
import os
import re
import sys

import numpy as np
import pandas as pd

from schrodinger.application.desmond import constants
from schrodinger.application.desmond.measurement import Measurement
# from schrodinger.application.scisol.packages.fep import fep_stats
from schrodinger.application.scisol.packages.fep import graph
from schrodinger.application.scisol.packages.fep_gui import unit_convert
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


# 1-letter standard amino acid residue code
INPUT_MUTATION_REGEX = re.compile(r'''
    ^
    (?P<chain>[A-Z])        # chain name
    [:-]                    # hyphen/dash or colon separator
    (?P<start>[A-Z])        # start residue 1-letter code
    (?P<resnum>\d+)         # numeric residue number
    (?P<inscode>[A-Za-z]?)  # optional insertion code
    (?:->?)?                # optional hyphen or mutations.txt-style arrow
    (?P<end>[A-Z]|\w{3})    # end residue 1- or 3-letter code (for NCAA)
    $
    ''', re.VERBOSE)


def _standardize_mutation_str(mut_str):
    '''Parse mutation from flexible input format to a tuple.'''
    if mut_str.upper() == 'WT':
        return 'WT'
    muts = mut_str.split(',')
    std_muts = []
    for m in muts:
        logger.debug(m)
        try:
            match = re.match(INPUT_MUTATION_REGEX, m).groupdict()
        except AttributeError:
            msg = f"Unable to parse mutation string: {mut_str}"
            raise ValueError(msg)
        std_mut = '{chain}:{start}{resnum}{inscode}->{end}'.format(**match)
        std_muts.append(std_mut.upper())
    return ','.join(std_muts)


def read_expt_file(inputfile):
    # colnames = ['mutations', 'KI', 'inequality']
    df = pd.read_csv(inputfile)
    logger.debug(df)
    l = df.values.tolist()
    d = {
        _standardize_mutation_str(x[0]): {
            'value': x[1],
            # 'inequality': None if pd.isna(x[2]) else x[2]
        }
        for x in l
    }
    logger.debug(d)
    return d


    # outdict = {}
    # with open(inputfile, 'r') as f:
    #     for line in f:
    #         # Skip empty lines and comments
    #         if len(line.strip()) == 0 or line[0] == "#":
    #             pass
    #         elif len(line.split(',')) > 0:
    #             key, value = [x.strip() for x in line.split(',')]
    #             try:
    #                 value = float(value)
    #             except ValueError:
    #                 # e.g. NA N/A NAN
    #                 value = None
    #             outdict[key] = value
    # return outdict


TITLE_REGEX = re.compile(r'''
    (?P<chain>[A-Za-z])-    # chain name + hyphen
    (?P<start>\w{3})        # start residue 3-letter code
    (?P<resnum>\d+)         # numeric residue number
    (?P<inscode>[A-Z]?)     # insertion code (is this always capitalized?)
    (?P<end>\w{3})          # end residue 3-letter code
    ''', re.VERBOSE)

AA_MAP = {
    'ALA': 'A',
    'CYS': 'C',
    'CYM': 'C',
    'ASP': 'D',
    'ASH': 'D',
    'GLU': 'E',
    'GLH': 'E',
    'PHE': 'F',
    'GLY': 'G',
    'HID': 'H',
    'HIE': 'H',
    'HIP': 'H',
    'ILE': 'I',
    'LYS': 'K',
    'LYN': 'K',
    'LEU': 'L',
    'MET': 'M',
    'ASN': 'N',
    'PRO': 'P',
    'GLN': 'Q',
    'ARG': 'R',
    'SER': 'S',
    'THR': 'T',
    'VAL': 'V',
    'TRP': 'W',
    'TYR': 'Y',
}


def convert_title_to_aa1(title):
    '''Convert from protein FEP node title format to single-letter AA format.'''
    sep = ','
    if title == "WT":
        return title
    muts3_in = title.split(sep)
    muts1_out = []
    for mut in muts3_in:
        m = re.match(TITLE_REGEX, mut).groupdict()

        # Handle noncanonical amino acids
        try:
            start = AA_MAP[m['start']]
        except KeyError:
            start = m['start']
        try:
            end = AA_MAP[m['end']]
        except KeyError:
            end = m['end']

        m1 = f"{m['chain']}:{start}{m['resnum']}{m['inscode']}->{end}"
        muts1_out.append(m1)
    title1 = sep.join(muts1_out)
    # print(title, '==>', title1)
    return title1


def main(fmp_in, data_file, units, relative=False):
    # input .fmp and data file
    # fmp_in = os.path.abspath(sys.argv[1])
    # data_file = os.path.abspath(sys.argv[2])
    fmp_basename, ext = os.path.splitext(fmp_in)
    fmp_out = fmp_basename + '_with_expt' + ext

    expt = read_expt_file(data_file)
    g = graph.Graph.deserialize(fmp_in)

    # all_stats = {}

    ref_exp_dg = None
    for n in g.nodes:
        title1 = convert_title_to_aa1(n.struc.title)
        logger.info(title1)

        val = None
        exp_dg = None

        try:
            val = expt[title1]['value']
        except KeyError:
            if title1 == "WT":
                print(f"WT not in experimental data, defaulting to ddG=0...")
                exp_dg = 0
            else:
                print(f"WARNING: {title1} not found in input data file..."
                      f"skipping node {n.struc.title}")
                continue

        # Convert to dG
        try:
            if val is not None:
                exp_dg = unit_convert.convert(
                    float(val), getattr(unit_convert.AffinityUnits, units)
                )
                assert not np.isnan(exp_dg)
        except ValueError:
            # Skip non-numeric values
            print(f"WARNING: could not convert value to float for {title1}: "
                  f"{val} ... skipping")
            continue
        except AssertionError:
            print(f"WARNING: skipping `nan` value for {title1}")
            continue

        try:
            binding_leg = n.get_leg_by_name(constants.PhysicalLegTypes.BINDING)
            binding_leg.exp_dg = Measurement(exp_dg)
            # folding_leg = n.get_leg_by_name(constants.PhysicalLegTypes.FOLDING)
            # folding_leg.exp_dg = Measurement(0)
        except AttributeError:
            # old graph or release, no physical legs
            n.exp_dg = Measurement(exp_dg)
        print(f"Added experimental value for node: {n.struc.title}")

        # TODO allow different reference node
        if n.struc.title == 'WT':
            ref_exp_dg = n.exp_dg
            print(f'  - Stored {n.struc.title} dG ({n.exp_dg}) '
                  'as reference value')

    if relative:
        print('\n→ Calculating relative dG values...')
        for n in g.nodes:
            if n.exp_dg is not None:
                # shift values without changing error magnitudes
                # print("Orig pred", n.pred_dg)
                # print("Orig expt", n.exp_dg)
                n.exp_dg = n.exp_dg - ref_exp_dg.val
                # n.pred_dg = n.pred_dg - ref_pred_dg.val
                # print("Final pred", n.pred_dg)
                # print("Final expt", n.exp_dg, '\n')


    # print('\n→ Calculating cycle closure...')
    g.calc_cycle_closure()

    # print('\n→ Calculating node stats...')
    # for n in g.nodes:
    #     if n.exp_dg is not None:
    #         # store stats
    #         node_stats = {
    #             'fep_error': n.pred_dg - n.exp_dg,
    #             'fep_abs_error': abs(n.pred_dg - n.exp_dg),
    #         }
    #         all_stats[n] = node_stats

    # report largest errors
    # largest_errors = sorted(
    #     all_stats.items(),
    #     key=lambda x: all_stats[x]['fep_abs_error'],
    #     reverse=True
    # )
    # print('  - Largest errors:\n')
    # [print(n, all_stats[n]['fep_abs_error']) for n, _ in largest_errors]

    if False:
        csv_out = f'{fmp_basename}_dG.csv'
        print(f'\n→ Writing {csv_out}...')
        with open(csv_out, 'w') as f:
            f.write('Ligand,exp_dg,pred_dg,pred_dg_error,label\n')
            for n in g.nodes:
                ligand = n.struc.title

                try:
                    exp_dg = n.exp_dg.val
                except AttributeError as e:
                    # NoneType
                    print(f'WARNING: no experimental dG for {ligand}')
                    exp_dg = ''

                try:
                    pred_dg = n.pred_dg.val
                    pred_dg_error = n.pred_dg.unc
                except AttributeError as e:
                    # None
                    print(f'WARNING: no predicted dG for {ligand}')
                    pred_dg = ''
                    pred_dg_error = ''

                label = n.short_id
                f.write(f'"{ligand}",{exp_dg},{pred_dg},{pred_dg_error},{label}\n')

    # print summary
    # print('\n→ Calculating summary statistics...\n')
    # results = fep_stats.calculate(g)
    # for item in results:
    #     print(item, results[item])

    print(f'\n→ Writing {fmp_out}...')
    g.write(fmp_out)

    # report largest errors


# TODO add `-m/--mutation-format` argument
# e.g. default `c:sri->e` == `{chain}:{start_aa}{resnum}{inscode}->{end_aa}`
# could use S/E to indicate 3-letter code
# could use R to indicate resnum + inscode
# TODO add `-s/--mutation-separator` argument (default ',')
# TODO add `-r/--reference-node` argument (default 'WT')
def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('fmp_in',
                        help='input .fmp file')
    parser.add_argument('data_file',
                        help='input file with experimental data')
    parser.add_argument('-u', '--units', default="KI_nM",
                        help='units of experimental data (one of: DG, KI_M, '
                             'KI_mM, KI_uM, KI_nM, PKI; default = KI_nM)')
    parser.add_argument('--relative', action='store_true',
                        help='use ddGs relative to WT/reference node')
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

    main(args.fmp_in, args.data_file, args.units, relative=args.relative)
