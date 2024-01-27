#!/bin/bash

# This script extracts pH vs. charge state #
# for the folded and unfolded states of    #
# coordinate file read by get_ss_charge.py #

# Steven Heaton
# 2024-01-27

source_dir="/project/target/complete"

find "$source_dir" -type f -name "*.log" -print0 | while IFS= read -r -d $'\0' file; do

    base_filename="$( basename -- "$file" )"
    protein_name="${base_filename%.*}"
    output_filename="${source_dir}/${protein_name}/apbs/${protein_name}_charge_values.tsv"

    awk '/Protein charge of folded and unfolded state as a function of pH/{flag=1; next} {if(flag && ++count<=142) {gsub(/[ \t]+/, "\t"); sub(/^\t/, ""); print}}' "$file" \
		> "$output_filename"

done
