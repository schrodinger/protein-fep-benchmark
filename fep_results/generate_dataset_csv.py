import argparse
import logging
from pathlib import Path
from pprint import pformat
import re

import numpy as np
import pandas as pd

from schrodinger.utils import log
from schrodinger.structure import RESIDUE_MAP_3_TO_1_LETTER


# Configure logging
DEFAULT_LOGGING_LEVEL = logging.WARNING
logger = log.get_output_logger(__name__)
logger.level = DEFAULT_LOGGING_LEVEL


RESULTS_PATH = Path('results')
RECALC_SUFFIX = '_recalc{}ns'
GROUP_DGS_SUFFIX = '_dGs_pH{}'


def separate_column(df: pd.DataFrame, col: str, into: list[str], sep: str,
                    remove: bool=True):
    # Split an input string column by `sep` into the columns listed in `into`.
    # For a non-string column, use its values for the first new column, and
    # set further columns to NaN.
    try:
        new_cols = df[col].str.split(pat=sep, expand=True)
        new_cols.columns = into
    except AttributeError:
        data_dict = {}
        for i in range(len(into)):
            if i == 0:
                data_dict[into[i]] = list(df[col])
            else:
                data_dict[into[i]] = [np.nan] * len(df[col])
        new_cols = pd.DataFrame.from_dict(data_dict, orient='columns')

    insert_index = df.columns.get_loc(col) + 1
    for i in range(len(new_cols.columns)):
        df.insert(insert_index + i, into[i], new_cols[into[i]])
    if remove:
        df = df.drop(col, axis=1)
    return df


def make_timepoint_df(system, time, dgs_csv):
    col_map = {
        'Node Group Title': 'node_group',
        'Node Title': 'mutation',
        'Node Experimental Binding dG': 'exp_dg_str',
        'Node Predicted Binding dG': 'naive_pred_dg_str',
        'dG Correction': 'dg_correction',
        'Node Group Predicted Binding dG': 'group_pred_dg_str',
    }
    dgs_df = (
        pd.read_csv(dgs_csv)
        .rename(columns=col_map)
    )
    # Add system and time columns
    dgs_df.insert(0, 'system', system)
    dgs_df.insert(1, 'time', system)

    dg_types = ['exp_dg', 'naive_pred_dg', 'group_pred_dg']
    for dg_type in dg_types:
        dgs_df = separate_column(dgs_df, f'{dg_type}_str',
                                 into=[dg_type, f'{dg_type}_unc'], sep=r"\+-")
    return dgs_df


def get_mutation_counts(mutations):
    mutation_counts = pd.Series([len(m) for m in mutations.str.split(",")])
    wt_indices = mutations[mutations=='WT'].index
    mutation_counts[wt_indices] = 0
    return mutation_counts


def generate_dataset(systems_csv):
    '''
    Generate the dataset.
    '''
    system_dict_list = pd.read_csv(systems_csv).to_dict(orient="records")

    # Read and merge all the system/timepoint dG csv files
    logger.info('\nReading individual system results'
                '\n---------------------------------')
    df = pd.DataFrame()
    for system_dict in system_dict_list: #[0:2]:
        system = system_dict['system']
        logger.info(f'\n{system}')
        logger.debug(pformat(system_dict, sort_dicts=False))
        logger.info(f'- reading dGs for {system_dict["orig_time"]} ns')
        orig_df = make_timepoint_df(system_dict['system'],
                                    system_dict['orig_time'],
                                    system_dict['orig_dgs_csv'])
        logger.info(f'- reading dGs for {system_dict["recalc_time"]} ns')
        recalc_df = make_timepoint_df(system_dict['system'],
                                      system_dict['recalc_time'],
                                      system_dict['recalc_dgs_csv'])
        logger.info(f'- merging into primary data frame')
        df = pd.concat([df, orig_df, recalc_df], axis=0, ignore_index=True)

    logger.info('\nProcessing merged data frame'
                '\n----------------------------')
    logger.debug(f'shape: {df.shape}')

    # Filter for single mutants
    logger.info('- Keeping only single mutants')
    df['n_mut'] = get_mutation_counts(df['mutation'])
    df = df[df['n_mut'] <= 1]
    logger.debug(f'shape: {df.shape}')

    # Discard missing experimental values
    logger.info('- Dropping missing experimental dGs')
    df = df.dropna(axis=0, subset=['exp_dg']).reset_index()
    logger.debug(f'shape: {df.shape}')


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('systems_csv', type=Path,
                        help='CSV file with systems metadata')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('--debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


def main(argv=None):
    '''Main workflow, optionally callable like subprocess with list of args.'''
    args = parse_args(argv)

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    generate_dataset(args.systems_csv)


if __name__ == '__main__':
    main()
