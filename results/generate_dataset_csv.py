import argparse
import logging
from pathlib import Path
from pprint import pformat, pprint
import re
from typing import Optional, Union

import numpy as np
import pandas as pd

from schrodinger.utils import log
from schrodinger.structure import RESIDUE_MAP_3_TO_1_LETTER
from schrodinger.application.scisol.fep import protein_fep_utils
from schrodinger.application.desmond.constants import WT_IDENTIFIER as WT


# Configure logging
DEFAULT_LOGGING_LEVEL = logging.WARNING
logger = log.get_output_logger(__name__)
logger.level = DEFAULT_LOGGING_LEVEL


# Filenames
RESULTS_PATH = Path('csv')
RECALC_SUFFIX = '_recalc{}ns'
GROUP_DGS_SUFFIX = '_dGs_pH{}'
FULL_DATASET_CSV = 'all_systems_full_dataset.csv'
BENCHMARK_CSV = 'benchmark_dataset.csv'
CASE_STUDIES_CSV = 'case_studies_dataset.csv'

# Data options
DG_PRECISION = 2
CASE_STUDY_SYSTEMS = ["3SKJ_HL", "5L6Y_HL", "6SBA_B"]

# Constants
NAIVE_RESCODES = ["ASP", "GLU", "HIE", "LYS"]
MUTATION_TYPE_LEVELS = ["WT", "neutral", "charged", "core-hopping"]
CHARGED_AA = [
    "ASP", "ASH",
    "GLU", "GLH",
    "HID", "HIE", "HIP", "HIS",
    "LYS", "LYN",
    "ARG"
]
CHARGED_AA1 = ["D", "E", "H", "K", "R"]
COREHOPPING_AA = ["PRO"]
COREHOPPING_AA1 = ["P"]
ALL_AA = RESIDUE_MAP_3_TO_1_LETTER.keys()

# Perturbation types
COREHOPPING = 'core-hopping'
CHARGED = 'charged'
NEUTRAL = 'neutral'

# ddG types
NAIVE = 'naive'
GROUP = 'group'

# System dict keys
ORIG_TIME = 'orig_time'
RECALC_TIME = 'recalc_time'
ORIG_DGS_CSV = 'orig_dgs_csv'
RECALC_DGS_CSV = 'recalc_dgs_csv'

# Column names
SYSTEM = 'system'
BASE_SYSTEM = 'base_system'
COMPONENT = 'component'
TIME = 'time'
NODE_GROUP = 'node_group'
MUTATION = 'mutation'
MUTATION1 = f'{MUTATION}1'
CHAIN = 'chain'
RESNUM = 'resnum'
INSCODE = 'inscode'
START_AA = 'start_aa'
END_AA = 'end_aa'
START_AA1 = f'{START_AA}1'
END_AA1 = f'{END_AA}1'
N_MUT = 'n_mut'
EXP_DG = 'exp_dg'
EXP_DG_UNC = 'exp_dg_unc'
DG_CORRECTION = 'dg_correction'
# DG_CORRECTION_UNC = 'dg_correction_unc'
NAIVE_PRED_DG = 'naive_pred_dg'
NAIVE_PRED_DG_UNC = 'naive_pred_dg_unc'
GROUP_PRED_DG = 'group_pred_dg'
GROUP_PRED_DG_UNC = 'group_pred_dg_unc'
WT_EXP_DG = 'wt_exp_dg'
WT_NAIVE_PRED_DG = 'wt_naive_pred_dg'
WT_GROUP_PRED_DG = 'wt_group_pred_dg'
NAIVE_PRED_DDG = 'naive_pred_ddg'
GROUP_PRED_DDG = 'group_pred_ddg'
EXP_DDG = 'exp_ddg'
DDG_TYPE = 'ddg_type'
PRED_DDG = 'pred_ddg'
ERR = 'err'
ABS_ERR = 'abs_err'
N_MICROSTATES = 'n_microstates'
TYPE = 'type'

# Used in matching WT and mutant dGs, for calculating ddG columns
WT_DG_JOIN_COLS = [SYSTEM, TIME]

# Column names used when matching protein_fep_utils.MUTATION_TITLE_RE
MUTATION_COLS = [CHAIN, START_AA, RESNUM, INSCODE, END_AA]

# Column names used when removing alternate protonation states
DEDUP_MERGE_COLS = [SYSTEM, MUTATION1, TIME]


def separate_column(df: pd.DataFrame, col: str, into: list[str],
                    pattern: Union[str, re.Pattern], remove: bool=True):
    # Split an input string column by `pattern` into the columns from `into`.
    # For a non-string column, use its values for the first new column, and
    # set further columns to NaN.
    try:
        new_cols = df[col].str.split(pat=pattern, expand=True)
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
        df = df.drop(columns=col)
    return df


def make_timepoint_df(system, time, dgs_csv):
    col_map = {
        'Node Group Title': NODE_GROUP,
        'Node Title': MUTATION,
        'Node Experimental Binding dG': f'{EXP_DG}_str',
        'Node Predicted Binding dG': f'{NAIVE_PRED_DG}_str',
        'dG Correction': DG_CORRECTION,
        'Node Group Predicted Binding dG': f'{GROUP_PRED_DG}_str',
    }
    dgs_df = (
        pd.read_csv(dgs_csv)
        .rename(columns=col_map)
    )
    # Add system and time columns
    dgs_df.insert(0, SYSTEM, system)
    dgs_df.insert(1, TIME, time)

    # Split ddG string columns
    # Protein FEP Groups script outputs Measurement strings as `{val}+-{unc}`
    dg_types = [EXP_DG, NAIVE_PRED_DG, GROUP_PRED_DG]
    for dg_col in dg_types:
        str_col = f'{dg_col}_str'
        unc_col = f'{dg_col}_unc'
        dgs_df = separate_column(dgs_df, str_col, into=[dg_col, unc_col],
                                 pattern=r"\+-", remove=True)
        dgs_df[dg_col] = dgs_df[dg_col].astype(float).round(DG_PRECISION)
        dgs_df[unc_col] = dgs_df[unc_col].astype(float).round(DG_PRECISION)

    return dgs_df


def get_mutation_counts(mutations):
    mutation_counts = pd.Series([len(m) for m in mutations.str.split(",")])
    wt_indices = mutations[mutations==WT].index
    mutation_counts[wt_indices] = 0
    return mutation_counts


def calculate_ddg_col(df: pd.DataFrame, dg_col: str,
                      ddg_col: Optional[str]=None,
                      join_by: Optional[list[str]]=None,
                      keep_wt: bool=False):
    if ddg_col is None:
        ddg_col = dg_col.replace('_dg', '_ddg')
    if join_by is None:
        join_by = WT_DG_JOIN_COLS
    wt_dg_col = f'wt_{dg_col}'
    wt_df = df[df[MUTATION] == WT]
    wt_df = wt_df.rename(columns={dg_col: wt_dg_col})
    wt_df = wt_df[[*WT_DG_JOIN_COLS, wt_dg_col]]
    out_df = df.merge(wt_df, how='left', left_on=join_by, right_on=join_by)
    out_df[ddg_col] = out_df[dg_col] - out_df[wt_dg_col]
    out_df[ddg_col] = out_df[ddg_col].round(2)
    if not keep_wt:
        out_df = out_df.drop(columns=wt_dg_col)
    return out_df


def deduplicate_protonation_states(df: pd.DataFrame,
                                   keep_rescodes: list[str]=NAIVE_RESCODES):
    '''
    Remove duplicate protonation states, keeping only rows with the specified
    residue codes in the END_AA column.
    '''
    n_microstates_df = (df.groupby(DEDUP_MERGE_COLS)
                          .agg({MUTATION: 'count'})
                          .rename(columns={MUTATION: N_MICROSTATES}))
    out_df = df.merge(n_microstates_df, how='left', left_on=DEDUP_MERGE_COLS,
                      right_on=DEDUP_MERGE_COLS)
    one_microstate = pd.Series(out_df[N_MICROSTATES] == 1)
    naive_rescode = pd.Series([(x in NAIVE_RESCODES) for x in out_df[END_AA]])
    keep = one_microstate | naive_rescode
    return out_df[keep]


def get_perturbation_type(start_aa, end_aa):
    '''
    Determine the perturbation typE based on the start and end residue codes.
    '''
    if (start_aa in COREHOPPING_AA) or (end_aa in COREHOPPING_AA):
        return 'core-hopping'
    elif (start_aa in CHARGED_AA) or (end_aa in CHARGED_AA):
        return 'charged'
    elif isinstance(start_aa, str) and isinstance(end_aa, str):
        return 'neutral'
    else:
        return 'WT'


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
        logger.info(f'\n{system_dict[SYSTEM]}')
        logger.debug(pformat(system_dict, sort_dicts=False))
        logger.info(f'- reading dGs for {system_dict[ORIG_TIME]} ns')
        orig_df = make_timepoint_df(system_dict[SYSTEM],
                                    system_dict[ORIG_TIME],
                                    system_dict[ORIG_DGS_CSV])
        logger.info(f'- reading dGs for {system_dict[RECALC_TIME]} ns')
        recalc_df = make_timepoint_df(system_dict[SYSTEM],
                                      system_dict[RECALC_TIME],
                                      system_dict[RECALC_DGS_CSV])
        logger.info(f'- merging into primary data frame')
        df = pd.concat([df, orig_df, recalc_df], axis=0, ignore_index=True)

    logger.info('\nProcessing merged data frame')
    logger.debug(f'    shape: {df.shape}')

    # Filter for single mutants
    logger.info('- Keeping only single mutants')
    df[N_MUT] = get_mutation_counts(df[MUTATION])
    df = df[df[N_MUT] <= 1]
    logger.debug(f'    shape: {df.shape}')

    # Discard missing experimental values
    logger.info('- Dropping missing experimental dGs')
    df = df.dropna(axis=0, subset=[EXP_DG]).reset_index(drop=True)
    logger.debug(f'    shape: {df.shape}')

    # Split mutations
    # Pandas Series.str.split() with regex produces an initial column for rows
    # that don't match the regex.  There is also a trailing column for some
    # reason which we need to remove.
    mutation_columns = ['_nomatch_', *MUTATION_COLS, '_extra_']
    df = separate_column(df, MUTATION, into=mutation_columns,
                         pattern=protein_fep_utils.MUTATION_TITLE_RE,
                         remove=False)
    df = df.drop(columns=['_nomatch_', '_extra_'])

    # Add 1-letter start/end residue codes and mutation1 column.values
    for restype in [START_AA, END_AA]:
        rescodes3 = df[restype].values.tolist()
        rescodes1 = [
            RESIDUE_MAP_3_TO_1_LETTER[res] if res is not None else None
            for res in rescodes3
        ]
        # Insert directly after 3-letter column.
        df.insert(df.columns.get_loc(restype) + 1, f'{restype}1', rescodes1)

    mutations1 = [protein_fep_utils.convert_node_title_3to1(mut)
                  for mut in df[MUTATION].values.tolist()]
    # Use colon chain separator
    mutations1 = [x.replace('-', ':') for x in mutations1]
    df.insert(df.columns.get_loc(MUTATION) + 1, MUTATION1, mutations1)

    # Add base_system and component
    df = separate_column(df, SYSTEM, into=[BASE_SYSTEM, COMPONENT],
                         pattern='_', remove=False)

    # Add ddG columns
    logger.info(f'- Calculating ddG columns (vs. {WT})')
    df = calculate_ddg_col(df, EXP_DG, keep_wt=True)
    df = calculate_ddg_col(df, NAIVE_PRED_DG, keep_wt=True)
    # df = calculate_ddg_col(df, DG_CORRECTION)
    df = calculate_ddg_col(df, GROUP_PRED_DG, keep_wt=True)

    # Remove WT rows
    logger.info(f'- Removing {WT} and self-mutation (titration) rows')
    df = df[df[MUTATION] != WT]
    # Remove titration mutations
    df = df[df[START_AA1] != df[END_AA1]]
    logger.debug(f'    shape: {df.shape}')

    logger.info('- Removing rows with redundant protonation states')
    df = deduplicate_protonation_states(df)
    logger.debug(f'    shape: {df.shape}')

    logger.info('- Making long w.r.t. pred. ddG type (group/naive)')
    # long_df = pd.melt(df, id_vars=[SYSTEM, MUTATION1, TIME],
    #                   value_vars=[NAIVE_PRED_DG, GROUP_PRED_DG])
    pivot_cols = [
        x for x in df.columns
        if x not in [NAIVE_PRED_DDG, GROUP_PRED_DDG]
    ]
    df = (df.rename(columns={NAIVE_PRED_DDG: NAIVE, GROUP_PRED_DDG: GROUP})
            .melt(id_vars=pivot_cols,
                  value_vars=[NAIVE, GROUP],
                  var_name=DDG_TYPE, value_name=PRED_DDG))
    logger.info('- Calculating err and abs_err')
    df[ERR] = round(df[PRED_DDG] - df[EXP_DDG], DG_PRECISION)
    df[ABS_ERR] = abs(df[ERR])


    # Add perturbation type
    logger.info('- Adding perturbation type')
    type_condition_list = [
        df[START_AA].isin(COREHOPPING_AA) | df[END_AA].isin(COREHOPPING_AA),
        df[START_AA].isin(CHARGED_AA) | df[END_AA].isin(CHARGED_AA),
        df[START_AA].isin(ALL_AA) & df[END_AA].isin(ALL_AA),
    ]
    type_choice_list = [COREHOPPING, CHARGED, NEUTRAL]
    df[TYPE] = np.select(type_condition_list, type_choice_list, default=np.nan)

    # Final sort
    logger.info('- Sorting sensibly')
    sort_columns = [SYSTEM, CHAIN, RESNUM, INSCODE, END_AA, TIME, DDG_TYPE]
    df = df.set_index(sort_columns).sort_index().reset_index()
    # Final column order
    final_columns = [
        SYSTEM,
        MUTATION,
        MUTATION1,
        CHAIN,
        START_AA,
        START_AA1,
        RESNUM,
        INSCODE,
        END_AA,
        END_AA1,
        BASE_SYSTEM,
        COMPONENT,
        TIME,
        NODE_GROUP,
        EXP_DG,
        EXP_DG_UNC,
        NAIVE_PRED_DG,
        NAIVE_PRED_DG_UNC,
        DG_CORRECTION,
        # DG_CORRECTION_UNC,
        GROUP_PRED_DG,
        GROUP_PRED_DG_UNC,
        N_MUT,
        WT_EXP_DG,
        WT_NAIVE_PRED_DG,
        WT_GROUP_PRED_DG,
        EXP_DDG,
        N_MICROSTATES,
        DDG_TYPE,
        PRED_DDG,
        ERR,
        ABS_ERR,
        TYPE,
    ]
    df = df[final_columns]
    # df.reindex(columns=final_columns)
    logger.debug(f'    shape: {df.shape}')

    # Write
    logger.info(f'\nWriting final output files:')

    df.to_csv(FULL_DATASET_CSV, index=False)
    logger.warning(f'- Wrote {FULL_DATASET_CSV}.')

    df[~df[SYSTEM].isin(CASE_STUDY_SYSTEMS)].to_csv(BENCHMARK_CSV, index=False)
    logger.warning(f'- Wrote {BENCHMARK_CSV}.')

    df[df[SYSTEM].isin(CASE_STUDY_SYSTEMS)].to_csv(CASE_STUDIES_CSV, index=False)
    logger.warning(f'- Wrote {CASE_STUDIES_CSV}.')

    # [print(f"'{x}',") for x in list(df.columns)]


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
