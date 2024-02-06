##################################################
# Config
##################################################
options(digits=4, scipen=5, width=80)
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggthemes)
    library(argparse)
    library(withr)  # avoid `could not find function "deferred_run"` error
})
# disable tidyverse deprecation warnings
rlang::local_options(lifecycle_verbosity = "quiet")


##################################################
# Constants / Lookups
##################################################

# gas constant
R <- 8.314 / 4184  # kcal mol-1 K-1

# Calculate energies
# We have to take some care with experimental inequalities.
ineq_lookup <- tribble(
    ~wt_ineq, ~mut_ineq, ~ddG_ineq,
    # If wt is exact, the mut inequality applies to ddG.
    "", "", "",
    "", "<", "<",
    "", ">", ">",
    # If mut is exact but wt is not, the ddG inequality is opposite that of wt.
    "<", "", ">",
    ">", "", "<",
    # If wt and mutant have the same inequality, we can conclude nothing.
    "<", "<", NULL,
    ">", ">", NULL,
    # If wt and mutant are opposite, the mut inequality applies to ddG.
    "<", ">", ">",
    ">", "<", "<"
)



##################################################
# Argument parsing
##################################################

# Parse command-line arguments
parse_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("systems", action="append", nargs="*",
                        help="PDB IDs for desired systems")
    parser$add_argument("-m", "--max-mutations", dest="max_mutations",
                        default=1, type="integer",
                        help="max mutations per case [default: %(default)s]")
    parser$add_argument("--min-cases", dest="min_cases",
                        default=10, type="integer",
                        help="min cases per system [default: %(default)s]")
    parser$add_argument("-i", "--input-file", dest="input_file",
                        default="skempi_v2.csv",
                        help="input CSV file [default: %(default)s]")
    parser$add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
                        help="don't write output files")
    parser$parse_args()
}



##################################################
# Functions
##################################################

# Convert SKEMPI single mutation to single-letter FEP node title syntax.
skempi_mut_to_fep_node_title <- function(x) {
    chain <- str_sub(x, 2, 2)
    start <- str_sub(x, 1, 1)
    num <- str_sub(x, 3, -2)
    end <- str_sub(x, -1, -1)
    str_c(chain, "-", start, num, end, sep="")
}


# Convert SKEMPI single mutation to single-letter FEP-style mutation syntax
# (with text arrow).
skempi_mut_to_fep_mutations1 <- function(x) {
    chain <- str_sub(x, 2, 2)
    start <- str_sub(x, 1, 1)
    num <- str_sub(x, 3, -2)
    end <- str_sub(x, -1, -1)
    str_c(chain, ":", start, num, "->", end, sep="")
}


# Convert SKEMPI single mutation to a named list.
skempi_mut_to_list <- function(x) {
    chain <- str_sub(x, 2, 2)
    start <- str_sub(x, 1, 1)
    num <- str_sub(x, 3, -2)
    end <- str_sub(x, -1, -1)
    list(chain=chain, start=start, resi=num, end=end)
}


# Split the SKEMPI pdb column into system and 2 sets of protein chain IDs
# Required columns: pdb
split_pdb_protein_chains <- function(df) {
    df %>%
        separate(
            pdb,
            c("system", "protein1_chains", "protein2_chains"),
            sep="_"
        )
}


# Convert from SKEMPI format to 1-letter FEP node title-like format
# Required columns: mutations
format_mutation_columns <- function(df) {
    df %>%
        mutate(
            # Temporary vector columns
            mutations_vec = str_split(mutations, ","),
            mutations_title_vec = lapply(mutations_vec, skempi_mut_to_fep_node_title),
            # mutations_fep_vec = lapply(mutations_vec, skempi_mut_to_fep_mutations1),

            # Convert to named list for future use
            mutations_list = lapply(mutations_vec, skempi_mut_to_list),
        ) %>%

        # Join vector columns rowwise
        rowwise() %>%
        mutate(
            n_mutations = length(mutations_vec),
            mutations_title = paste(mutations_title_vec, collapse=","),
            # mutations_fep = paste(mutations_fep_vec, collapse=","),
        ) %>%
        ungroup() %>%
        select(-mutations_vec, -mutations_title_vec)
}


# Restrict to SKEMPI-categorized locations
# Required column: location
limit_to_locations <- function(df, locations){
    df %>%
        rowwise() %>%
        mutate(
            location_vec = str_split(location, ","),
            location_uniq_vec = list(sort(unique(location_vec))),
            location_uniq = paste(location_uniq_vec, collapse=","),
            has_only_desired_locations = all(location_vec %in% locations),
        ) %>%
        ungroup() %>%
        filter(has_only_desired_locations) %>%
        select(-location_vec, -location_uniq_vec, -has_only_desired_locations)
}


# Remove cases with mutations on both sides of the interface
# Required columns: mutations_list, protein1_chains, protein2_chains
exclude_cross_interface_mult_muts <- function(df) {
    df %>%
        rowwise() %>%
        mutate(
            mut_chains = unique(c(mutations_list["chain"])),
            mutates_protein1 = !all(is.na(str_match(protein1_chains, mut_chains))),
            mutates_protein2 = !all(is.na(str_match(protein2_chains, mut_chains))),
            mutates_both_proteins = mutates_protein1 & mutates_protein2,
        ) %>%
        ungroup() %>%
        filter(!mutates_both_proteins) %>%
        select(-mutates_protein1, -mutates_protein2, -mutates_both_proteins) %>%
        rowwise() %>%
        mutate(
            all_mut_chains = paste(unique(mut_chains), collapse=","),
        ) %>%
        ungroup()
}


# Only keep cases measured using specific methods
# Required column: method
limit_to_methods <- function(df, methods){
    filter(df, method %in% methods)
}


# Required columns: temp, affinity_{wt,mut}_parsed
calculate_dG_and_ddG <- function(df) {
    df %>%
    mutate(
        wt_dG = R * temp * log(affinity_wt_parsed),
        mut_dG = R * temp * log(affinity_mut_parsed),
        expt_ddG = mut_dG - wt_dG,
    )
}


# Average multiple measurements
# Required columns: system, protein1_chains, protein2_chains mutations_title,
#                   n_mutations, location_uniq, expt_ddG, method
calculate_mean_ddG <- function(df) {
    df %>%
    group_by(
        system,
        protein1_chains, protein2_chains,
        mut_chains,
        # all_mut_chains,
        mutations_title,
        n_mutations,
        location_uniq,
    ) %>%
    summarize(
        .groups = "keep",
        expt_ddG_mean = mean(expt_ddG),
        expt_ddG_sd = sd(expt_ddG),
        N_meas = n(),
        methods = str_c(method, sep=","),
    ) %>%
    ungroup() %>%
    # select(-mutations_list) %>%
    distinct() %>%
    replace_na(list(expt_ddG_sd = 0))
}


# Write mutations for a system to one or two CSV files
# Required columns:
#   system, expt_ddG, expt_ddG_sd, protein1_chains, protein2_chains,
#   mut_in_p1_only, mut_in_p2_only, mutations_title, N_meas, location_uniq
extract_system_to_csv <- function(df, pdbid, precision=2) {
    xdf <- df %>%
        filter(system == pdbid) %>%
        mutate(
            expt_ddG = round(expt_ddG, digits=precision),
            expt_ddG_sd = round(expt_ddG_sd, digits=precision),
        )

    # Get the chain IDs for each protein component
    p1 <- sort(unique(xdf$protein1_chains))
    p2 <- sort(unique(xdf$protein2_chains))

    # Set output filenames
    p1_fn <- paste(pdbid, p1, "expt.csv", sep="_")
    p2_fn <- paste(pdbid, p2, "expt.csv", sep="_")

    # First protein component
    p1_df <- xdf %>% filter(mut_in_p1_only == T)
    if (nrow(p1_df) > 0) {
        p1_df %>%
            select(
                variant=mutations_title,
                expt_ddG,
                # expt_ddG_sd,
                # N_meas,
                # location_uniq
            ) %>%
            distinct() %>%
            write_csv(p1_fn)
        cat(paste("\nWrote", p1_fn, "with", nrow(p1_df), "rows."))
    }

    # 2nd component
    p2_df <- xdf %>% filter(mut_in_p2_only == T)
    if (nrow(p2_df) > 0) {
        p2_df %>%
            select(
                variant=mutations_title,
                expt_ddG,
                # expt_ddG_sd,
                # N_meas,
                # location_uniq
            ) %>%
            distinct() %>%
            write_csv(p2_fn)
        cat(paste("\nWrote", p2_fn, "with", nrow(p2_df), "rows."))
    }
}



##################################################
# Main Workflow
##################################################

# Parse args
args <- parse_args()

# Load SKEMPI v2 data
raw_df <- read_delim(args$input_file, delim=";",
                     col_types=cols(.default = col_character()))

df <- raw_df %>%
    # Keep only needed columns and rename
    select(
        pdb = `#Pdb`,
        location = `iMutation_Location(s)`,
        mutations = `Mutation(s)_PDB`,
        affinity_mut_M = `Affinity_mut (M)`,
        affinity_mut_parsed = Affinity_mut_parsed,
        affinity_wt_M = `Affinity_wt (M)`,
        affinity_wt_parsed = Affinity_wt_parsed,
        temp = Temperature,
        method = Method,
    ) %>%

    # Create columns for
    split_pdb_protein_chains() %>%
    format_mutation_columns() %>%

    # Remove cases involving protein core and non-interface surface
    limit_to_locations(c("COR", "SUP", "RIM")) %>%

    # Exclude cases with mutations on both sides of the interface, which are
    # not supported by FEP.
    exclude_cross_interface_mult_muts() %>%

    # Use only SPR and ITC data
    limit_to_methods(c("SPR", "ITC")) %>%

    # Convert values to numeric
    mutate(
        # Binding affinities
        affinity_mut_parsed = as.numeric(affinity_mut_parsed),
        affinity_wt_parsed = as.numeric(affinity_wt_parsed),

        # Temperature (treat assumed as given)
        temp = as.numeric(str_replace(temp, "\\(assumed\\)", "")),
    ) %>%

    # Extract inequality operators from the affinity value strings
    mutate(
        wt_ineq = if_else(str_sub(affinity_wt_M, 1, 1) %in% c("<", ">"),
                          str_sub(affinity_wt_M, 1, 1),
                          ""),
        mut_ineq = if_else(str_sub(affinity_mut_M, 1, 1) %in% c("<", ">"),
                           str_sub(affinity_mut_M, 1, 1),
                           ""),
    ) %>%

    # Merge in the lookup to create ddG_ineq column
    left_join(ineq_lookup, by=c("wt_ineq", "mut_ineq")) %>%

    # Flatten lists (not sure why these are lists!)
    # mutate(
    #     ddG_ineq = str_c(ddG_ineq),
    # ) %>%

    # For now, remove inexact measurements.
    filter(ddG_ineq == "") %>%

    # Calculate dG and ddG values, averaging for cases with >1 measurement.
    calculate_dG_and_ddG() %>%
    calculate_mean_ddG() %>%

    # Discard data points with large experimental standard deviations.
    # filter(expt_ddG_sd <= 1)

    # Discard cases with too many mutations
    filter(n_mutations <= args$max_mutations) %>%

    # Select columns for our final "clean" experimental data frame.
    select(
        system,
        protein1_chains,
        protein2_chains,
        mut_chains,
        # all_mut_chains,
        mutations_title,
        n_mutations,
        location_uniq,
        # ddG_ineq,
        expt_ddG=expt_ddG_mean,
        expt_ddG_sd,
        N_meas,
        methods
    ) %>%

    # Drop incomplete records.
    drop_na() %>%

    # Track which component has each mutated position
    # TODO: Move this earlier
    rowwise() %>%
    mutate(
        mut_in_p1 = any(str_detect(protein1_chains, unlist(mut_chains))),
        mut_in_p2 = any(str_detect(protein2_chains, unlist(mut_chains))),
        mut_in_p1_only = mut_in_p1 & !mut_in_p2,
        mut_in_p2_only = mut_in_p2 & !mut_in_p1,
        mut_in_both = mut_in_p1 & mut_in_p2,
    ) %>%
    ungroup()


# Maybe eventually output some plots
# summary(df)
# qplot(methods, data=df)
# qplot(expt_ddG, data=df)
# qplot(expt_ddG_sd, data=df)
# qplot(n_mutations, data=df)
# qplot(N_meas, data=df)
# qplot(location_uniq, data=df)



##################################################
# Summarize Systems
##################################################

# Other potential filters
# min_dynamic_range <- 2.5
# interesting_ddG_cutoff <- 1

# Limit to requested systems if requested
if (length(unlist(args$systems)) > 0) {
    systems_df <- df %>%
        filter(system %in% args$systems)
} else {
    # Otherwise include all systems
    systems_df <- df
}

# Calculate a summary of the systems satisfying the selection criteria
systems_df <- systems_df %>%
    group_by(system, protein1_chains, protein2_chains) %>%
    summarize(
        .groups = "keep",
        N = n(),
        # dynamic_range = max(expt_ddG) - min(expt_ddG),
        # min_ddG = min(expt_ddG),
    ) %>%
    filter(
        N >= args$min_cases
        # & dynamic_range >= min_dynamic_range
        # & min_ddG <= interesting_ddG_cutoff
    ) %>%
    arrange(-N)

cat("\nSKEMPI Systems Summary")
cat("\n======================")
knitr::kable(systems_df, format="simple")

cat(paste("\nTotal number of systems:", nrow(systems_df)))
cat(paste("\nTotal number of cases:", sum(systems_df$N)))



##################################################
# Write output files
##################################################
if (args$dry_run) {
    cat("\n\nDry run - not writing output files.\n")
} else {
    cat("\n\nWriting output files...\n")
    for (system in unique(systems_df %>% pull(system))) {
        extract_system_to_csv(df, system)
    }
}

cat("\n\nDone.\n\n")
