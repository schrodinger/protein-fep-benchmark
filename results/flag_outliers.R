#!/usr/bin/env Rscript

# Global options
options(digits=4, scipen=5, width=80)

# Load libraries quietly
suppressWarnings(
    suppressMessages({
        library(tidyverse)
        library(argparse)
    })
)



##################################################
# Setup
##################################################
options(digits=4, scipen=5, width=80)
suppressPackageStartupMessages({
    library(tidyverse)
    library(argparse)
})



##################################################
# Defaults
##################################################

# Tunable parameters
DEFAULT_FSASA_CUTOFF <- 0.1
DEFAULT_CHARGED_NEIGHBORS_CUTOFF <- 6
DEFAULT_ABS_DELTA_SIZE_CUTOFF <- 5
DEFAULT_INTERACTION_FREQ_CUTOFF <- 0.5
DEFAULT_DDG_CUTOFF_FEP <- 2.5

# DEFAULT_GLY_NEIGHBOR_CUTOFF <- 30
DEFAULT_OUT_CSV <- "flagged_cases.csv"
DEFAULT_CHARGED_FLAG_CORR_TERM_OUTFILE <- "charged_flag_correction_term.txt"

# Meta parameters
FSASA_SCALE <- 20
PARTIAL_FSASA_SCALE <- 10
ABS_DELTA_SIZE_SCALE <- 1
DDG_SCALE_FEP <- 1
CHARGED_NEIGHBORS_SCALE <- 1
INTERACTION_FREQ_SCALE <- 10


# Score cutoffs (generally 0.5 ** <number of continuous parameters>)
CUTOFF_SBB_XINT <- 0.5 ** 1
CUTOFF_SBB_SOL_BURIED <- 0.5 ** 1
CUTOFF_INTRODUCES_BURIED_CHARGE <- 0.5 ** 3
CUTOFF_CHARGE_FLIP <- 0.5 ** 3
CUTOFF_BREAKS_PI_INTERACTION <- 0.5 ** 3
CUTOFF_BURIED_AROMATIC <- 0.5 ** 2
CUTOFF_PROLINE <- 0.5
CUTOFF_BURIED_GLYCINE <- 0.75  # higher cutoff for stricter burial requirements



##################################################
# Constants
##################################################

AROMATIC_1LETTER <- c('H', 'F', 'Y', 'W')



##################################################
# Argument parsing
##################################################

# Parse command-line arguments
parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("features_csv",
                        type="character",
                        help="input features CSV file")
    parser$add_argument("--fsasa-cutoff",
                        type="numeric",
                        dest="fsasa_cutoff",
                        default=DEFAULT_FSASA_CUTOFF,
                        help="fractional SASA cutoff [default: %(default)s]")
    parser$add_argument("--ddg-cutoff",
                        type="numeric",
                        dest="ddg_cutoff_fep",
                        default=DEFAULT_DDG_CUTOFF_FEP,
                        help="high-energy ddG cutoff [default: %(default)s]")
    parser$add_argument("--abs-delta-size-cutoff",
                        type="numeric",
                        dest="abs_delta_size_cutoff",
                        default=DEFAULT_ABS_DELTA_SIZE_CUTOFF,
                        help="absolute delta size cutoff (# of non-H atoms) [default: %(default)s]")
    parser$add_argument("--charged-neighbors-cutoff",
                        type="numeric",
                        dest="charged_neighbors_cutoff",
                        default=DEFAULT_CHARGED_NEIGHBORS_CUTOFF,
                        help="charged neighbors cutoff (# residues within 5 Å) [default: %(default)s]")
    parser$add_argument("--interaction-freq-cutoff",
                        type="numeric",
                        dest="interaction_freq_cutoff",
                        default=DEFAULT_INTERACTION_FREQ_CUTOFF,
                        help="trajectory-based interaction frequency cutoff [default: %(default)s]")
    # parser$add_argument("--gly-neighbor-cutoff", dest="gly_neighbor_cutoff",
    #                     default=DEFAULT_GLY_NEIGHBOR_CUTOFF,
    #                     help="glycine neighbor cutoff (# of non-H atoms within 5 Å) [default: %(default)s]")
    parser$add_argument("--charged-flag-corr", dest="charged_flag_corr", type="numeric", default=0,
                        help="empirical charged flag correction term [default: %(default)s]")
    parser$add_argument(
        "-o", "--out-csv",
        dest="out_csv",
        default=DEFAULT_OUT_CSV,
        help="output CSV filename [default: %(default)s]"
    )
    parser$add_argument("--debug", dest="debug", action="store_true",
                        help="output debugging info")
    parser$parse_args()
}



##################################################
# Functions
##################################################

# Return the score for x > x0 from a continuous sigmoid curve.
#
# Parameters for the logistic growth curve:
#   x = value to score
#   x0 = threshold value
#   k = scale_factor (controls steepness of curve at midpoint)
score_greater_than <- function(x, x0, k=1) {
    1 / (1 + exp(-k * (x - x0)))
}


# Return the score for x < x0 from a continuous sigmoid curve.
#
# Same parameters as `score_greater_than`.
score_less_than <- function(x, x0, k=1) {
    1 - score_greater_than(x, x0, k=k)
}



##################################################
# Main Workflow
##################################################

# Parse args
args <- parse_args()

# Default or user-defined params
FSASA_CUTOFF <- args$fsasa_cutoff
DDG_CUTOFF_FEP <- args$ddg_cutoff_fep
ABS_DELTA_SIZE_CUTOFF <- args$abs_delta_size_cutoff
CHARGED_NEIGHBORS_CUTOFF <- args$charged_neighbors_cutoff
INTERACTION_FREQ_CUTOFF <- args$interaction_freq_cutoff
# GLY_NEIGHBOR_CUTOFF <- args$gly_neighbor_cutoff
OUT_CSV = args$out_csv
CHARGED_FLAG_CORR_TERM_OUTFILE <- DEFAULT_CHARGED_FLAG_CORR_TERM_OUTFILE

# Derived params
DDG_CUTOFF_PRIME <- 10 * DDG_CUTOFF_FEP
PARTIAL_FSASA_CUTOFF <- 2.5 * FSASA_CUTOFF



# Load df
features_df <- read_csv(args$features_csv, show_col_types=F)


# Initial processing of features_df
df <- features_df %>%
    # result features
    mutate(
        outlier_1.0 = abs_err > 1.0,
        outlier_1.5 = abs_err > 1.5,
        outlier_2.0 = abs_err > 2.0,
    ) %>%

    # rounding (eventually move this to the merge script)
    mutate(
        err = round(err, 2),
        abs_err = round(abs_err, 2),
    ) %>%

    # derived features
    mutate(
        # SASA
        cpx_buried = score_less_than(start_cpx_fsasa, FSASA_CUTOFF,
                                     FSASA_SCALE),
        sol_buried = score_less_than(start_sol_fsasa, FSASA_CUTOFF,
                                     FSASA_SCALE),
        mut_cpx_buried = score_less_than(end_cpx_fsasa, FSASA_CUTOFF,
                                         FSASA_SCALE),
        cpx_mostly_buried = score_less_than(start_cpx_fsasa,
                                            PARTIAL_FSASA_CUTOFF,
                                            PARTIAL_FSASA_SCALE),
        sol_mostly_buried = score_less_than(start_sol_fsasa,
                                            PARTIAL_FSASA_CUTOFF,
                                            PARTIAL_FSASA_SCALE),
        mut_cpx_mostly_buried = score_less_than(end_cpx_fsasa,
                                                PARTIAL_FSASA_CUTOFF,
                                                PARTIAL_FSASA_SCALE),
        # start_small_delta_fsasa = score_less_than(
        #     abs(start_delta_fsasa) / start_sol_fsasa,
        #     0.5,
        #     1
        # ),
        # end_small_delta_fsasa = score_less_than(
        #     abs(end_delta_fsasa) / end_sol_fsasa,
        #     0.5,
        #     1
        # ),

        # Charged environment
        highly_charged_env = score_greater_than(charged_neighbors_5.0,
                                                CHARGED_NEIGHBORS_CUTOFF,
                                                CHARGED_NEIGHBORS_SCALE),

        # ddG
        destabilizing_pred = score_greater_than(pred_ddg, DDG_CUTOFF_FEP,
                                                DDG_SCALE_FEP),

        # residue type
        wt_aromatic = as.integer(start_aa1 %in% AROMATIC_1LETTER),
        mut_aromatic = as.integer(end_aa1 %in% AROMATIC_1LETTER),
        proline = as.integer(str_detect(mutation, "PRO")),
        wt_glycine = as.integer(start_aa1 == "G"),
        glycine_loop = wt_glycine & res_ss_str == "loop",
        # gly_crowded = as.integer(wt_glycine) * score_greater_than(
        #         start_n_cpx_ha_neighbors, GLY_NEIGHBOR_CUTOFF),

        # residue size
        large_size_change = score_greater_than(abs_delta_size,
                                               ABS_DELTA_SIZE_CUTOFF,
                                               ABS_DELTA_SIZE_SCALE),

        # charge
        introduces_charge = as.integer(start_formal_charge == 0 &
                                           end_formal_charge != 0),
        flips_charge = as.integer(sign(start_formal_charge) *
                                      sign(end_formal_charge)) < 0,

        # pi interactions
        wt_pipi = score_greater_than(max_sum_wt_pipi_freq,
                                     INTERACTION_FREQ_CUTOFF,
                                     INTERACTION_FREQ_SCALE),
        mut_pipi = score_greater_than(max_sum_mut_pipi_freq,
                                      INTERACTION_FREQ_CUTOFF,
                                      INTERACTION_FREQ_SCALE),
        wt_picat = score_greater_than(max_sum_wt_picat_freq,
                                      INTERACTION_FREQ_CUTOFF,
                                      INTERACTION_FREQ_SCALE),
        mut_picat = score_greater_than(max_sum_mut_picat_freq,
                                       INTERACTION_FREQ_CUTOFF,
                                       INTERACTION_FREQ_SCALE),
        wt_xint_pi_hbond = as.integer(start_n_cpx_pi_hbonds >
                                          start_n_sol_pi_hbonds),
        mut_xint_pi_hbond = as.integer(end_n_cpx_pi_hbonds >
                                           end_n_sol_pi_hbonds),
    ) %>%

    # Limit to results at the longest available simulation time  # FIXME
    filter(time == max(time, na.rm=T)) %>%


    ### FLAGS ###

    mutate(
        # Breaks cross-interface salt bridge at buried site in complex
        flag_sbb_xint = coalesce(
            as.integer(breaks_cross_int_salt_bridge) *
                pmax(cpx_buried, highly_charged_env),
            0),


        # Breaks salt bridge at buried site in the solvent leg
        flag_sbb_sol_buried = coalesce(
            as.integer(breaks_salt_bridge & !breaks_cross_int_salt_bridge)
                * sol_buried,
            0),


        # Introduces buried charge at buried site & destabilizing
        flag_introduces_buried_charge = coalesce(
            introduces_charge * cpx_mostly_buried * mut_cpx_buried *
                destabilizing_pred,
            0),


        # Flips charge at mostly buried site
        flag_charge_flip = coalesce(
            flips_charge * cpx_mostly_buried * mut_cpx_mostly_buried *
                destabilizing_pred,
            0),


        # Breaks pipi or picat interaction
        breaks_pipi = wt_pipi * wt_aromatic *
            pmax((1 - mut_pipi), (1 - mut_aromatic)),
        breaks_picat = wt_picat *
            as.integer(wt_aromatic == 1 | start_formal_charge == 1) *
            pmax(
                (1 - mut_picat),
                as.integer(wt_aromatic & !mut_aromatic),
                as.integer(start_formal_charge == 1 &
                               end_formal_charge < 1)
            ),
        breaks_pi_hbond = as.integer(wt_xint_pi_hbond & !mut_xint_pi_hbond),
        flag_breaks_pi_interaction = coalesce(
            pmax(breaks_pipi, breaks_picat) * cpx_buried,
            0),


        # Mutates aromatic residue at buried site with large size change
        flag_buried_aromatic = coalesce(
            wt_aromatic * cpx_buried * large_size_change *
                destabilizing_pred,
            0),


        # Mutation to or from proline
        flag_proline = coalesce(as.integer(proline), 0),


        # Mutation from glycine at buried site
        flag_buried_glycine = coalesce(
            as.integer(wt_glycine) * as.integer(!glycine_loop) * cpx_buried,
            0),

        # flag_non_contact = coalesce(
        #     start_non_contact * end_non_contact
        # )
    ) %>%


    # Normalize flags
    mutate(
        norm_flag_sbb_xint = flag_sbb_xint / CUTOFF_SBB_XINT,
        norm_flag_sbb_sol_buried = flag_sbb_sol_buried / CUTOFF_SBB_SOL_BURIED,
        # norm_flag_breaks_pi_interaction = flag_breaks_pi_interaction /
        #     CUTOFF_BREAKS_PI_INTERACTION,
        norm_flag_buried_aromatic = flag_buried_aromatic /
            CUTOFF_BURIED_AROMATIC,
        norm_flag_proline = flag_proline / CUTOFF_PROLINE,
        norm_flag_buried_glycine = flag_buried_glycine / CUTOFF_BURIED_GLYCINE,
        norm_flag_introduces_buried_charge = flag_introduces_buried_charge /
            CUTOFF_INTRODUCES_BURIED_CHARGE,
        norm_flag_charge_flip = flag_charge_flip / CUTOFF_CHARGE_FLIP,
    ) %>%


    # Collect all flags
    mutate(
        # if any normalized flag > 1 the case is flagged
        flagged = pmax(
            norm_flag_sbb_xint,
            norm_flag_sbb_sol_buried,
            # norm_flag_breaks_pi_interaction,
            norm_flag_buried_aromatic,
            norm_flag_proline,
            # norm_flag_buried_glycine,
            norm_flag_introduces_buried_charge,
            norm_flag_charge_flip,
            na.rm=T
        ) > 1,
    )


##################################################
# Determine empirical correction
##################################################

# Combine all `flag_*` column information into a single output column
separate_flags <- function(X) {
    X %>%
        pivot_longer(
            cols = matches("norm_flag_.*"),
            names_to = "flag_name",
            names_prefix = "norm_flag_",
            values_to = "flag_value",
        )
}

charged_flags_df <- tribble(
    ~flag_name, ~charged_flag, ~affected_leg,
    "sbb_xint", TRUE, "complex",
    "sbb_sol_buried", TRUE, "solvent",
    "introduces_buried_charge", TRUE, "complex",
    "charge_flip", TRUE, "complex",
    "breaks_pi_interaction", FALSE, "complex",
    "buried_aromatic", FALSE, "complex",
    "proline", FALSE, "both",
    "buried_glycine", FALSE, "complex",
)

correction_df <- df %>%
    filter(flagged) %>%
    separate_flags() %>%
    filter(flag_value > 1) %>%
    select(system, mutation1, exp_ddg, pred_ddg, err, abs_err,
           flag_name, flag_value) %>%
    left_join(charged_flags_df, by="flag_name") %>%
    mutate(leg_err = if_else(affected_leg == "solvent", -err, err)) %>%
    filter(charged_flag==TRUE) %>%
    group_by(charged_flag) %>%
    summarize(
        n = n(),
        mean_signed_leg_err = mean(leg_err),
        sd_signed_leg_err = sd(leg_err)
    )

complex_charged_flag_names <- charged_flags_df %>%
    filter(charged_flag & affected_leg == "complex") %>%
    pull(flag_name) %>%
    unlist()
solvent_charged_flag_names <- charged_flags_df %>%
    filter(charged_flag & affected_leg == "solvent") %>%
    pull(flag_name) %>%
    unlist()
complex_charged_flag_names_to_detect <- paste(complex_charged_flag_names, collapse="|")
solvent_charged_flag_names_to_detect <- paste(solvent_charged_flag_names, collapse="|")

if (args$charged_flag_corr == 0) {
    charged_flag_correction_term <- round(correction_df$mean_signed_leg_err, 1)
    charged_flag_correction_term_sd <- round(correction_df$sd_signed_leg_err, 1)
    cat("Calculated charged flag correction term: ", charged_flag_correction_term,
        "+-", charged_flag_correction_term_sd, "\n")
} else {
    charged_flag_correction_term <- args$charged_flag_corr
    cat("Using given charged flag correction term: ", charged_flag_correction_term, "\n")
}
complex_charged_flag_correction <- -1 * charged_flag_correction_term
solvent_charged_flag_correction <- 1 * charged_flag_correction_term

# Write the value to the correction file
corr_file_handle <- file(CHARGED_FLAG_CORR_TERM_OUTFILE)
writeLines(paste0(charged_flag_correction_term), corr_file_handle)
close(corr_file_handle)
cat("→ Wrote charged flag correction term (", charged_flag_correction_term,
    ") to file: ", CHARGED_FLAG_CORR_TERM_OUTFILE, ".\n", sep="")



##################################################
# Prepare output table
##################################################

# Combine all `flag_*` column information into a single output column
consolidate_flags <- function(X) {
    X %>%
        pivot_longer(
            cols = matches("^norm_flag_.*"),
            names_to = "flag_name",
            names_prefix = "norm_flag_",
            values_to = "flag_value",
        ) %>%
        filter(flag_value > 1) %>%
        group_by(across(-matches("flag_.*|norm_flag_.*"))) %>%
        summarize(
            .groups = "keep",
            flags = paste(flag_name, collapse="|"),
        ) %>%
        ungroup()
}


# Debug: write full flag features
if (args$debug) {
    full_flags_csv <- paste0('full_dev_features_and_flags.csv')
    df %>%
        select(
            system,
            mutation1,
            start_cpx_fsasa,
            cpx_buried,
            cpx_mostly_buried,
            matches("ddg$"),
            matches("err$"),
            matches("^flag"),
        ) %>%
        write_csv(full_flags_csv)
}


# Report and write flagged cases CSV.
flagged_df <- df %>%
    filter(flagged) %>%
    consolidate_flags() %>%
    mutate(
       charged_flag_correction = case_when(
               str_detect(flags, complex_charged_flag_names_to_detect) ~ complex_charged_flag_correction,
               str_detect(flags, solvent_charged_flag_names_to_detect) ~ solvent_charged_flag_correction,
               TRUE ~ 0),
    ) %>%
    select(system, mutation, mutation1, matches("^flag_.*"), flags, charged_flag_correction)

cat(paste0("Flagged ", nrow(flagged_df), " cases out of ", nrow(df), " total cases.\n"))
write_csv(flagged_df, OUT_CSV)





##################################################
# Stats
##################################################

# Get the full path of this script file
# https://stackoverflow.com/a/45389384/501277
thisFile <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
        # 'source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
}


stats_out_csv <- "param_fitting_data.csv"
script_file <- thisFile()  # rstudioapi::getSourceEditorContext()$path
# script_hash <- digest::digest(script_file, file=T)
stats_df <- as_tibble_row(
        lst(
            script_hash = digest::digest(script_file, file=T),
            timestamp_ms = round(as.numeric(Sys.time())*1000),
            input_csv = args$features_csv,
            fsasa_cutoff=FSASA_CUTOFF,
            ddg_cutoff = DDG_CUTOFF_FEP,
            abs_delta_size_cutoff = ABS_DELTA_SIZE_CUTOFF,
            charged_neighbors_cutoff = CHARGED_NEIGHBORS_CUTOFF,
            interaction_freq_cutoff = INTERACTION_FREQ_CUTOFF,
            fsasa_scale = FSASA_SCALE,
            partial_fsasa_scale = PARTIAL_FSASA_SCALE,
            ddg_scale = DDG_SCALE_FEP,
            charged_neighbors_scale = CHARGED_NEIGHBORS_SCALE,
            interaction_freq_scale = INTERACTION_FREQ_SCALE,
            n_cases = nrow(df),
            n_outliers_2.0 = nrow(filter(df, outlier_2.0)),
            n_outliers_1.5 = nrow(filter(df, outlier_1.5)),
            n_outliers_1.0 = nrow(filter(df, outlier_1.0)),
            n_flagged_outliers = nrow(filter(df, flagged & outlier_2.0)),
            n_flagged_grey_area = nrow(filter(df, flagged & outlier_1.0 & !outlier_2.0)),
            n_flagged_good_pred = nrow(filter(df, flagged & !outlier_1.0)),
            n_unflagged_outliers = nrow(filter(df, !flagged & outlier_2.0)),
            n_unflagged_grey_area = nrow(filter(df, !flagged & outlier_1.0 & !outlier_2.0)),
            n_unflagged_good_pred = nrow(filter(df, !flagged & !outlier_1.0)),
        )
    )
write_csv(stats_df, stats_out_csv, append=file.exists(stats_out_csv))
cat("→ Wrote stats to file: ", stats_out_csv, ".\n", sep="")



cat(paste0("→ Wrote output flags file: ", OUT_CSV, ".\nDone.\n\n"))
