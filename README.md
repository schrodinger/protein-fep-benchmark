# Protein Binding FEP+ Benchmark

Data and analysis scripts for Sampson, J. M. et al. "Robust prediction of relative binding energies for protein-protein complex mutations using free energy perturbation calculations."  J. Mol. Biol. 436, 168640 (2024). [https://doi.org/10.1016/j.jmb.2024.168640]()


The following analysis was performed on the systems in the paper.

1. Prepare and process the experimental data, downloading the SKEMPI2 data file if necessary.

    ```bash
    ( cd experimental_data && ./prepare_data.sh )
    ```

2. Place original 100 ns `out.fmp` files into `results/fmp`.  Due to size limitations, these files are not included in this repository.

3. Edit the `results/systems.csv` file to include the relevant metadata for each system.

4. Compile the FEP results dataset, including experimental data and derived data columns.

    ```bash
    ( cd results && ./compile_dataset.sh -v )
    ```

5. Flag probable outliers from the .fmp files.

    ```bash
    export SCHRODINGER=/path/to/Schrodinger_Suite_2024-2>
    (
        cd results

        # Calculate the .fmp features (this may take a while!)
        time ./calculate_features.sh -v

        # Merge the features and results into a single dataset for stats
        Rscript merge_features_and_results.R

        # Do the same for testing the case studies dataset
        Rscript merge_features_and_results.R --case-studies

        # Flag probable outliers in the benchmark dataset
        # (writes `charged_flag_correction_term.txt`)
        Rscript flag_outliers.R \
            benchmark_merged_features_results.csv \
            -o benchmark_flagged_cases.csv

        # Flag probable outliers in the case study dataset using
        # `charged_flag_correction_term.txt` from the benchmark dataset
        Rscript flag_outliers.R \
            case_studies_merged_features_results.csv \
            -o case_studies_flagged_cases.csv
    )
    ```
