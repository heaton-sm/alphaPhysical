#!/bin/bash

# This script extracts the stability and   #
# pH vs. charge state for the folded and   #
# unfolded states of a coordinate file     #
# read by get_ss_charge_stability.py       #

# Steven Heaton
# 2024-02-01

source_dir="/project/target/complete"

find "$source_dir" -type f -name "*.log" -print0 | while IFS= read -r -d $'\0' propka_out; do

    base_filename="$( basename -- "${propka_out}" )"
    protein_name="${base_filename%.*}"
    charge_filename="${source_dir}/${protein_name}/apbs/${protein_name}_charge_values.tsv"
	stability_filename="${source_dir}/${protein_name}/apbs/${protein_name}_stability_values.tsv"

    awk '/Protein charge of folded and unfolded state as a function of pH/{flag=1; next} {if(flag && ++count<=142) {gsub(/[ \t]+/, "\t"); sub(/^\t/, ""); print}}' "${propka_out}" \
		> "${charge_filename}"
    awk '/Free energy of   folding/{flag=1; next} {if(flag && ++count<=16) {gsub(/[ \t]+/, "\t"); sub(/^\t/, ""); print}}' "${propka_out}" \
		> "${stability_filename}"
	sed -i '1s/^/pH\tfree_energy\n/' "${stability_filename}"

done
