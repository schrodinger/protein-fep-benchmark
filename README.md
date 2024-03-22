# Protein Binding FEP+ Benchmark
Data and scripts for the Protein Binding FEP+ benchmark manuscript (Sampson, et al., in preparation).


1. Prepare and process the experimental data, downloading the SKEMPI2 data file if necessary.
```bash
( cd experimental_data && ./prepare_data.sh )
```

2. Compile the FEP results dataset, including experimental data and derived data columns.
```bash
( cd fep_results && ./compile_dataset.sh )
```
