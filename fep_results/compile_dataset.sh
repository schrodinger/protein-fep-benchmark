#!/usr/bin/env bash

# Inputs
reprocess_time=10

# Paths
results_dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
process_script="$results_dir/process_systems.py"
generate_dataset_script="$results_dir/generate_dataset_csv.py"
systems_csv="$results_dir/systems.csv"
systems_out_csv="$results_dir/systems_out.csv"

# Logging
now="$(date +%Y-%m-%d-%H-%M-%S)"
log_dir="${results_dir}/log"
mkdir -p "$log_dir"
log="${log_dir}/${now}_compile_dataset.log"

main () {
    (

    # Make sure we're running from this directory
    cd "$results_dir"

    # Process the out.fmp files
    "$SCHRODINGER/run" "$process_script" \
        -i "$systems_csv" \
        -t "$reprocess_time" \
        -v "$@" \
        | tee "$log"

    # Compile the individual fmp results into a single dataset
    $SCHRODINGER/run "$generate_dataset_script" \
        "$systems_out_csv" -v \
        | tee 

    )
}

main "$@"
