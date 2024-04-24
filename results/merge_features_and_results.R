
options(digits=4, scipen=5, width=80)

# Global options
options(digits=4, scipen=5, width=80)

# Load libraries quietly
suppressWarnings(
    suppressMessages({
        library(tidyverse)
        library(ggthemes)
        library(argparse)
        # library(ggbeeswarm)
        # library(caret)
    })
)

##################################################
# Arguments
##################################################

# Parse command-line arguments
parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--case-studies", action="store_true",
                        help="use case studies dataset")
    # parser$add_argument("-i", "--input-file", required=TRUE,
    #                     help="input file")
    # parser$add_argument("-v", "--verbose", action="store_true",
    #                     help="increase verbosity")
    # parser$add_argument("--debug", action="store_true",
    #                     help="print debugging info")
    parser$parse_args()
}


# Parse args
args <- parse_args()


##################################################
#
##################################################

# Input files
if (args$case_studies) {
    fep_results_csv_file <- "case_studies_dataset.csv"
    out_csv <- "case_studies_merged_features_results.csv"
} else {
    fep_results_csv_file <- "benchmark_dataset.csv"
    out_csv <- "benchmark_merged_features_results.csv"
}
# mutpred_merged_out_csv_file <- "../csv/mutpred_merged_out.csv"
features_csv_files <- Sys.glob("features/features*.csv")



##################################################
# Functions
##################################################

arr_str_to_numeric <- function(arr_str) {
    extracted <- str_extract_all(arr_str, "\\d+\\.\\d+")[[1]]
    as.numeric(extracted)
}

get_arr_index_nearest_val <- function(arr, val) {
    # index of the value nearest zero for the array shifted by -val
    which.min(abs(arr - val))
}


interpolate_from_arr_str <- function(x_str, y_str, xout) {
    x <- arr_str_to_numeric(x_str)
    y <- arr_str_to_numeric(y_str)
    approx(x, y, xout)$y
    # idx <- get_arr_index_nearest_val(time_arr, time)
    # val_arr[idx]
}



##################################################
# Main Workflow
##################################################

# Load and process FEP results
fep_raw <- read_csv(fep_results_csv_file, show_col_types=F)
fep_df <- fep_raw %>%
    filter(ddg_type == "group") %>%
    select(
        # Prune unneeded columns
        -n_mut, -n_microstates, -matches("^naive"), -ddg_type,
    ) %>%
    mutate(
        time = factor(time, levels=c(10, 100)),
    ) %>%
    mutate(
        # Derived values
        err = pred_ddg - exp_ddg,
        abs_err = abs(err),
        outlier_2 = abs_err > 2,
        outlier_1.5 = abs_err > 1.5,
    ) %>%
    select(
        system, mutation, mutation1, type, time, exp_ddg, pred_ddg, err, abs_err
    )

aromatic_aa1 <- c("F", "Y", "W")
nonpolar_aa1 <- c("A", "V", "I", "L", "M")
polar_aa1 <- c("C", "S", "T", "N", "Q")
charged_aa1 <- c("D", "E", "H", "K", "R")
glypro_aa1 <- c("G", "P")

aa1_restype_map <- tibble(
    aa1 = c(aromatic_aa1, nonpolar_aa1, polar_aa1, charged_aa1, glypro_aa1),
    restype = c(
        rep("aromatic", length(aromatic_aa1)),
        rep("nonpolar", length(nonpolar_aa1)),
        rep("polar", length(polar_aa1)),
        rep("charged", length(charged_aa1)),
        rep("glypro", length(glypro_aa1))
    )
)


# Load mutpred results
# mutpred_df <- read_csv(mutpred_merged_out_csv_file, show_col_types=F) %>%
#     separate_wider_regex(
#         Mutations,
#         patterns=c(
#             chain="[A-Z]", ":", position="\\d+[A-Z]?",
#             "\\(", start="[A-Z]{3}", "->", end="[A-Z]{3}", "\\)"
#         )
#     ) %>%
#     mutate(
#         mutation = str_c(chain, "-", start, position, end, sep="")
#     ) %>%
#     relocate(system, mutation) %>%
#     select(-chain, -start, -position, -end)

# mutpred_df
# features_raw %>% filter(n0_aa1 == "WT") %>% select(system, mutation) %>% nrow()
# mutpred_df %>% select(system, mutation) %>% nrow()

# Load and process .fmp-derived features
features_raw <- do.call(bind_rows, lapply(features_csv_files, read_csv, show_col_types=F))
features_df <- features_raw %>%
    # Only use direct mutations (not cycle closure mutations to same endpoint)
    filter(n0_aa1 == "WT") %>%
    # Derived edge features
    mutate(
        mutation1 = n1_aa1,
        delta_size = end_size - start_size,
        abs_delta_size = abs(delta_size),
        proline = start_aa1=="P" | end_aa1=="P",
        glycine = start_aa1=="G" | end_aa1=="G",
        # time_arr = str_arr_to_numeric(time_arr),
        # complex_dg_arr = str_arr_to_numeric(complex_dg_arr),
        # solvent_dg_arr = str_arr_to_numeric(solvent_dg_arr),
        # ddg_arr = str_arr_to_numeric(ddg_arr),

        # Changes in pi-pi and pi-cation interaction frequencies
        delta_max_pipi_freq = max_mut_pipi_freq - max_wt_pipi_freq,
        delta_sum_pipi_freq = sum_mut_pipi_freq - sum_wt_pipi_freq,
        delta_max_picat_freq = max_mut_picat_freq - max_wt_picat_freq,
        delta_sum_picat_freq = sum_mut_picat_freq - sum_wt_picat_freq,

    ) %>%

    # Add restype columns
    left_join(aa1_restype_map, by=c("start_aa1"="aa1")) %>%
    rename(start_restype=restype) %>%
    left_join(aa1_restype_map, by=c("end_aa1"="aa1")) %>%
    rename(end_restype=restype) %>%

    # Get 1ns energy
    rowwise() %>%
    mutate(
        ddg_1ns = interpolate_from_arr_str(time_arr, ddg_arr, 1),
    ) %>%
    ungroup() %>%

    # left_join(mutpred_df, by=join_by(system, mutation)) %>%

    # Consolidate protonation microstates
    group_by(system, mutation1) %>%
    summarize(
        .groups = "keep",

        start_aa1 = first(start_aa1),
        end_aa1 = first(end_aa1),

        start_restype = first(start_restype),
        end_restype = first(end_restype),

        res_ss = first(res_ss),
        res_ss_str = first(res_ss_str),
        is_charged = any(is_charged),

        # Only one of the microstates should be charged, but normalize with
        # `sign` just in case
        start_formal_charge = sign(sum(start_formal_charge)),
        end_formal_charge = sign(sum(end_formal_charge)),
        delta_formal_charge = end_formal_charge - start_formal_charge,
        abs_delta_formal_charge = abs(delta_formal_charge),

        # if any of the edges involving the mutation break a salt
        # bridge, mark it as breaking a salt bridge
        breaks_salt_bridge = any(breaks_salt_bridge),
        breaks_cross_int_salt_bridge = any(breaks_cross_int_salt_bridge),

        # Keep size info (same for all microstates)
        start_size = first(start_size),
        end_size = first(end_size),
        delta_size = first(delta_size),
        abs_delta_size = first(abs_delta_size),

        # Use the most-buried fSASA values
        start_cpx_fsasa = min(start_cpx_fsasa),
        end_cpx_fsasa = min(end_cpx_fsasa),
        start_sol_fsasa = min(start_sol_fsasa),
        end_sol_fsasa = min(end_sol_fsasa),
        start_cpx_sc_fsasa = min(start_cpx_sc_fsasa),
        end_cpx_sc_fsasa = min(end_cpx_sc_fsasa),
        start_sol_sc_fsasa = min(start_sol_sc_fsasa),
        end_sol_sc_fsasa = min(end_sol_sc_fsasa),

        # Change in fSASA upon binding for start/end residues
        start_delta_fsasa = start_cpx_fsasa - start_sol_fsasa,
        end_delta_fsasa = end_cpx_fsasa - end_sol_fsasa,
        start_delta_sc_fsasa = start_cpx_sc_fsasa - start_sol_sc_fsasa,
        end_delta_sc_fsasa = end_cpx_sc_fsasa - end_sol_sc_fsasa,

        # energy change
        max_delta_prime_energy = max(delta_prime_energy),

        # max_delta_affinity = max(`delta Affinity`),
        # min_delta_affinity = min(`delta Affinity`),
        # max_delta_stability = max(`delta Stability`),
        # min_delta_stability = min(`delta Stability`),

        # FEP energy
        max_ddg_1ns = max(ddg_1ns),

        # pipi frequency
        max_wt_pipi_freq = max(max_wt_pipi_freq),
        max_mut_pipi_freq = max(max_mut_pipi_freq),
        max_sum_wt_pipi_freq = max(sum_wt_pipi_freq),
        max_sum_mut_pipi_freq = max(sum_mut_pipi_freq),

        # picat frequency
        max_wt_picat_freq = max(max_wt_picat_freq),
        max_mut_picat_freq = max(max_mut_picat_freq),
        max_sum_wt_picat_freq = max(sum_wt_picat_freq),
        max_sum_mut_picat_freq = max(sum_mut_picat_freq),

        # Delta max pipi frequency
        max_delta_max_pipi_freq = max(delta_max_pipi_freq),
        min_delta_max_pipi_freq = min(delta_max_pipi_freq),
        largest_delta_max_pipi_freq = if_else(
            max_delta_max_pipi_freq > abs(min_delta_max_pipi_freq),
            max_delta_max_pipi_freq,
            min_delta_max_pipi_freq,
        ),
        abs_largest_delta_max_pipi_freq = abs(largest_delta_max_pipi_freq),

        # Delta sum pipi frequency
        max_delta_sum_pipi_freq = max(delta_sum_pipi_freq),
        min_delta_sum_pipi_freq = min(delta_sum_pipi_freq),
        largest_delta_sum_pipi_freq = if_else(
            max_delta_sum_pipi_freq > abs(min_delta_sum_pipi_freq),
            max_delta_sum_pipi_freq,
            min_delta_sum_pipi_freq
        ),
        abs_largest_delta_sum_pipi_freq = abs(largest_delta_sum_pipi_freq),

        # Delta max picat frequency
        max_delta_max_picat_freq = max(delta_max_picat_freq),
        min_delta_max_picat_freq = min(delta_max_picat_freq),
        largest_delta_max_picat_freq = if_else(
            max_delta_max_picat_freq > abs(min_delta_max_picat_freq),
            max_delta_max_picat_freq,
            min_delta_max_picat_freq
        ),
        abs_largest_delta_max_picat_freq = abs(largest_delta_max_picat_freq),

        # Delta sum picat frequency
        max_delta_sum_picat_freq = max(delta_sum_picat_freq),
        min_delta_sum_picat_freq = min(delta_sum_picat_freq),
        largest_delta_sum_picat_freq = if_else(
            max_delta_sum_picat_freq > abs(min_delta_sum_picat_freq),
            max_delta_sum_picat_freq,
            min_delta_sum_picat_freq
        ),
        abs_largest_delta_sum_picat_freq = abs(largest_delta_sum_picat_freq),

        # Change in # pi Hbonds
        start_n_cpx_pi_hbonds = max(start_n_cpx_pi_hbonds),
        start_n_sol_pi_hbonds = max(start_n_sol_pi_hbonds),
        end_n_cpx_pi_hbonds = max(end_n_cpx_pi_hbonds),
        end_n_sol_pi_hbonds = max(end_n_sol_pi_hbonds),

        # Use the max charged neighbors at each cutoff
        charged_neighbors_5.0 = max(charged_neighbors_5.0),
        charged_neighbors_7.5 = max(charged_neighbors_7.5),
        charged_neighbors_10.0 = max(charged_neighbors_10.0),

        # Use the max heavy atom neighbors for each state
        start_n_cpx_ha_neighbors = max(start_n_cpx_ha_neighbors),
        end_n_cpx_ha_neighbors = max(end_n_cpx_ha_neighbors),
        start_n_sol_ha_neighbors = max(start_n_sol_ha_neighbors),
        end_n_sol_ha_neighbors = max(end_n_sol_ha_neighbors),
    )


# Merge the results and the features
df <- fep_df %>%
    left_join(features_df, by=join_by(system, mutation1)) %>%
    # Handle missing direct edges (e.g. charged <-> proline)
    mutate(
        breaks_salt_bridge = if_else(is.na(breaks_salt_bridge), FALSE, breaks_salt_bridge)
    ) %>%
    relocate(system, exp_ddg, mutation1) %>%
    write_csv(out_csv)
