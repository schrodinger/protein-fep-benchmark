##################################################
# Config
##################################################
options(digits=4, scipen=5, width=80)
suppressPackageStartupMessages({
    library(tidyverse)
    library(argparse)
    library(knitr)
    library(readxl)
    library(withr)  # avoid `could not find function "deferred_run"` error
})
# disable tidyverse deprecation warnings
rlang::local_options(lifecycle_verbosity = "quiet")


##################################################
# Constants / Lookups
##################################################

PRECISION <- 2
WT <- "WT"

# Gas constant
R <- 8.314 / 4184  # kcal mol-1 K-1

# Room temperature
temp <- 298.15



##################################################
# Argument parsing
##################################################

# Parse command-line arguments
parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("-i", "--input-file", dest="input_file",
                        required=T,
                        help="input Excel .xlsx file")
    parser$add_argument("systems", action="append", nargs="*",
                        help="PDB IDs for desired systems (tabs in input file)")
    parser$add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
                        help="don't write output files")
    parser$parse_args()
}



##################################################
# Main Workflow
##################################################

# Parse args
args <- parse_args()

cat(paste0("\nReading data from ", args$input_file, "\n"))

if (args$dry_run) {
    cat("\n\nDry run - not writing output files.\n")
} else {
    cat("\nWriting output files...\n")
}

for (pdbid in args$systems) {
    pdbid_df <- read_excel(args$input_file, sheet=pdbid) %>%
        mutate(
            system = str_c(pdbid, "_", mutation_chains),
        ) %>%
        filter(measurement == "Kd") %>%
        group_by(pdbid, system, variant) %>%
        summarize(
            .groups = "keep",
            affinity_mut = exp(mean(log(value))),  # geometric mean
        ) %>%
        ungroup() %>%
        left_join(
            by="pdbid",
            filter(., variant == "WT") %>%
                group_by(pdbid) %>%
                summarize(affinity_wt = mean(affinity_mut)) %>%
                ungroup() %>%
                select(pdbid, affinity_wt)
        ) %>%
        mutate(
            expt_ddG = R * temp * log(affinity_mut / affinity_wt),
            expt_ddG = round(expt_ddG, PRECISION),
        ) %>%
        filter(variant != WT) %>%
        select(system, variant, expt_ddG) %>%
        arrange(system, variant) %>%
        unique()


    for (sys in unique(pdbid_df$system)) {
        out_csv <- paste0(sys, "_expt.csv")
        out_df <- filter(pdbid_df, system == sys | variant == WT) %>%
            select(variant, expt_ddG)
        if (args$dry_run) {
            print(knitr::kable(out_df, format="simple"))
            cat("\nDRY RUN - no file was written.")
        } else {
            write_csv(out_df, out_csv)
            cat(paste("\nWrote", out_csv, "with", nrow(out_df), "rows."))
        }
    }
}

cat("\n\nDone\n\n")

