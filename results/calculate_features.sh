#!/usr/bin/env bash

fmp_dir="fmp"
fmps=$(cd "$fmp_dir" && find * -name "*with_expt.fmp" -maxdepth 0)

# echo $fmps | tr -s ' ' '\n' | wc -l

# Extract features from .fmp files
for fmp in $fmps; do
	fmp_path="${fmp_dir}/${fmp}"
	echo "############################################################"
	echo "# FMP: $fmp"
	# Pass any command line args (like -v or --debug) to the script
	$SCHRODINGER/run get_fmp_edge_features.py "$fmp_path" "$@"
	echo
done
