import argparse
import logging
import os
from pathlib import Path
from pprint import pformat
import re
import subprocess
import sys

import pandas as pd
from colorama import Fore

from schrodinger.application.scisol.fep import graph
from schrodinger.utils import log


# Configure logging
DEFAULT_LOGGING_LEVEL = logging.WARNING
logger = log.get_output_logger(__name__)
logger.level = DEFAULT_LOGGING_LEVEL


DEFAULT_MODEL_DG_PATTERN = 'model_dG_{0}.csv'
DEFAULT_OUT_DIR = Path("csv")

ADD_EXPT_SCRIPT = 'add_expt_to_fmp.py'
REANALYZE_SCRIPT = 'reanalyze_fmp.py'
PROTEIN_GROUPS_SCRIPT ='protein_fep_groups.py'


def parse_systems_csv(csv):
    '''Parse systems info csv and return a list of dicts.'''
    df = pd.read_csv(csv)
    return df.to_dict('records')


def all_outfiles_exist(*outfiles):
    '''
    Return True if all expected outfiles exist, False otherwise.
    '''
    for outfile in outfiles:
        if not outfile.is_file():
            return False
    return True


def skip_step(msg):
    '''
    Log a message to indicate the current step is being skipped.
    '''
    logger.info(f'{msg}...{Fore.YELLOW}skipping{Fore.RESET}')


def srun(args):
    '''Run a command with $SCHRODINGER/run via subprocess.'''
    schrod = os.environ['SCHRODINGER']
    subprocess.run([f'{schrod}/run'] + args, check=True)


def args_append_debug_verbosity(args):
    '''
    Return modified args with --debug or -v flags appended if appropriate.
    '''
    if logger.level <= logging.DEBUG:
        args.append('--debug')
    elif logger.level <= logging.INFO:
        args.append('-v')
    return args


def run_add_expt(fmp, expt, units='DG', force=False):
    '''Add experimental data to an .fmp and return the path to the new .fmp'''
    # script info
    script = ADD_EXPT_SCRIPT
    info_msg = f'adding experimental data from {expt} to {fmp}...'

    # output filename
    fmp_suffix = 'with_expt'
    fmp_base, _ = os.path.splitext(fmp)
    fmp_new = Path(f'{fmp_base}_{fmp_suffix}.fmp')

    if all_outfiles_exist(fmp_new) and not force:
        skip_step(f'Output file <{fmp_new}> already exists')
        return fmp_new

    # assemble command
    args = [
        Path(script).expanduser(),
        fmp,
        expt,
        '-u', units,
        '--relative',
    ]
    args = args_append_debug_verbosity(args)

    srun(args)
    assert(all_outfiles_exist(fmp_new))
    return fmp_new


def run_reanalyze_fmp(fmp, t, force=False):
    '''
    Recalculate dG values for the .fmp file at time `t`.

    Returns the path to the new .fmp file.
    '''
    # script info
    script = REANALYZE_SCRIPT

    # output filename
    fmp_suffix = f'recalc{t}ns'
    fmp_base, _ = os.path.splitext(fmp)
    fmp_new = Path(f'{fmp_base}_{fmp_suffix}.fmp')

    if all_outfiles_exist(fmp_new) and not force:
        skip_step(f'Output file <{fmp_new}> already exists')
        return fmp_new

    # assemble command
    args = [
        script,
        fmp,
        '-t', str(t),
        '-o', fmp_new,
    ]
    args = args_append_debug_verbosity(args)

    srun(args)
    assert(all_outfiles_exist(fmp_new))
    return fmp_new


def run_fep_groups(fmp, model_dg_csv, ph, out_dir=DEFAULT_OUT_DIR, force=False):
    '''
    Apply protein FEP groups correction to the .fmp file at the given pH.

    Returns the paths to the new dGs and pKas .csv file.

    '''
    fmp_base, _ = os.path.splitext(os.path.basename(fmp))
    out_dgs_csv = out_dir / f'{fmp_base}_dGs_pH{ph}.csv'
    out_pkas_csv = out_dir / f'{fmp_base}_pKas_pH{ph}.csv'

    outfiles = [out_dgs_csv, out_pkas_csv]
    if all_outfiles_exist(*outfiles) and not force:
        skip_step('Output files already exist')
        return outfiles

    args = [
        PROTEIN_GROUPS_SCRIPT,
        fmp,
        '-model_csv', model_dg_csv,
        '-ph', str(ph),
    ]
    srun(args)

    # Protein_fep_groups.py currently always writes to current working dir, so
    # we have to move the files to the result directory.
    for outfile in outfiles:
        os.rename(os.path.basename(outfile), outfile)

    assert(all_outfiles_exist(*outfiles))
    return outfiles


def detect_release(fmp):
    g = graph.Graph.deserialize(fmp)
    # Assume release used for first edge complex leg was used for entire graph
    first_sid_text = str(list(g.edges)[0].get_leg_by_name('complex'))
    # Match string like "ReleaseVersion=2021-2"
    pattern = r'''\d{4}-\d'''
    release = re.search(pattern, first_sid_text).group()
    # release = release_text.split("=")[1]
    logger.debug(f"{fmp} was run with release = {release}")
    return release


def write_output_systems_csv(result_dicts, fn="systems_out.csv"):
    '''
    Write the results to the output csv.

    '''
    for r in result_dicts:
        for k, v in r.items():
            if isinstance(v, Path):
                r[k] = v.as_posix()
    logger.debug('result_dicts')
    logger.debug(pformat(result_dicts, sort_dicts=False))
    df = pd.DataFrame(result_dicts)
    df.to_csv(fn, index=False)
    logger.info(f'Wrote {fn}')
    return fn


def process_system(d, recalc_time, out_dir=DEFAULT_OUT_DIR, force=False):
    '''Calculate the results for the system.'''
    # Release and model_dgs file
    release = detect_release(d['orig_fmp'])
    assert release == d['release']  # sanity check
    model_dg_csv = DEFAULT_MODEL_DG_PATTERN.format(release)
    logger.debug(f'using model_dG file = {model_dg_csv}')

    # Add expt data
    # logger.info(f'add experimental data')
    fmp_with_expt = run_add_expt(d['orig_fmp'], d['expt_csv'], force=force)

    # Reanalyze at new time
    logger.info(f'reanalyze fmp at {recalc_time}ns')
    recalc_fmp = run_reanalyze_fmp(fmp_with_expt, recalc_time, force=force)
    # recalc_fmp = run_reanalyze_fmp(d['orig_fmp'], recalc_time, force=force)

    # Groups correction
    logger.info(f'FEP groups calculation on original .fmp')
    orig_dgs_csv, orig_pkas_csv = run_fep_groups(
        fmp_with_expt, model_dg_csv, d['expt_ph'], out_dir=out_dir,
        #d['orig_fmp'], model_dg_csv, d['expt_ph'], out_dir=out_dir,
        force=force
    )

    logger.info(f'FEP groups calculation on recalculated .fmp')
    recalc_dgs_csv, recalc_pkas_csv = run_fep_groups(
        recalc_fmp, model_dg_csv, d['expt_ph'], out_dir=out_dir, force=force,
    )
    return {
        **d,
        'model_dg_csv': model_dg_csv,
        'fmp_with_expt': fmp_with_expt,
        'recalc_time': recalc_time,
        'recalc_fmp': recalc_fmp,
        'orig_dgs_csv': orig_dgs_csv,
        'orig_pkas_csv': orig_pkas_csv,
        'recalc_dgs_csv': recalc_dgs_csv,
        'recalc_pkas_csv': recalc_pkas_csv,
    }


def main(systems_csv, recalc_time, out_dir=DEFAULT_OUT_DIR, force=False):
    '''Main workflow'''
    systems_dicts = parse_systems_csv(systems_csv)

    result_dicts = []
    for d in systems_dicts:  # [:1]  # slice here for testing
        logger.info(f'\nWorking on {d["system"]}'
                    f'\n-----------{len(str(d["system"]))*"-"}')
        result = process_system(d, recalc_time, out_dir=out_dir, force=force)
        result_dicts.append(result)

    logger.info(f'\nFinishing up'
                f'\n------------')

    systems_base, csv_ext = os.path.splitext(systems_csv)
    out_systems_csv = f'{systems_base}_out{csv_ext}'
    write_output_systems_csv(result_dicts, fn=out_systems_csv)

    logger.info("Done.")


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-csv',
                        type=Path,
                        default='systems.csv',
                        help='CSV file with info for systems to be processed')
    parser.add_argument('-t', '--recalc-time',
                        type=int,
                        default=10,
                        help='reanalyze_fmp timepoint (ns)')
    parser.add_argument('-o', '--out-dir',
                        type=Path,
                        default=DEFAULT_OUT_DIR,
                        help='output directory for writing result files')
    parser.add_argument('-f', '--force',
                        action='store_true',
                        help='overwrite existing output files')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='increase verbosity')
    parser.add_argument('--debug',
                        action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    main(args.input_csv, args.recalc_time, out_dir=args.out_dir,
         force=args.force)
