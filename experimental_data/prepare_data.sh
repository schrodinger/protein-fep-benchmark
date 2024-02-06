#!/usr/bin/env bash

print_divider() {
    echo -e "\n#############################################################"
}

print_section_header() {
    print_divider
    echo -e "\n→ $1"
}


##################################################
# SKEMPI Systems
##################################################
skempi_csv="skempi_v2.csv"
skempi_url="https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv"

# Download the SKEMPI database CSV if necessary
if [[ ! -f "$skempi_csv" ]]; then
    print_section_header "Downloading SKEMPI CSV file..."
    curl -o "$skempi_csv" "$skempi_url"
    echo -e "\e  - Saved as $skempi_csv\n"
fi

# Extract the curated data points for the selected systems
print_section_header "Preparing SKEMPI system experimental CSV files"
Rscript prepare_skempi_systems.R \
    -i "$skempi_csv" \
    1DVF 1IAR 1JRH 1KTZ 1OGA 1VFB


print_divider


