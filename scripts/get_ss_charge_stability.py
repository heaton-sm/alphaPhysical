#!/usr/bin/env python3

# This script loads .pdb sequences and     #
# computes the percentage of alpha helix   #
# and beta sheet content using Pymol, then #
# calculates electrostatics using APBS.    #

# Steven Heaton
# 2024-02-01

import os
import subprocess
import pdb2pqr
import numpy as np
import concurrent.futures
from pymol import cmd
from natsort import natsorted

# function for reading coordinate file and extracting chain length and secondary structure content
def calculate_secondary_structure(pdb_file):
    cmd.load(pdb_file)
    atoms = cmd.get_model("n. ca").atom
    secstruct = [i.ss for i in atoms]

    if len(secstruct) == 0:
        helix_content = 0.0
        sheet_content = 0.0
    else:
        helix_content = 100.0 * secstruct.count("H") / len(secstruct)
        sheet_content = 100.0 * secstruct.count("S") / len(secstruct)

    chain_length = len(set(atom.resi for atom in atoms))
    cmd.reinitialize()
    return helix_content, sheet_content, chain_length

# function for running pdb2pqr
def run_pdb2pqr(pdb_file, output_pqr):
    command = ["pdb2pqr", "--ff=PARSE", "--whitespace", "--titration-state-method=propka", "--with-ph=7", pdb_file, output_pqr]
    subprocess.run(command)

# function for writing .pqr file using pdb2pqr for APBS run
def generate_apbs_input(output_pqr, output_apbs_in):
    with open(output_apbs_in, "w") as f:
        f.write("read\n")
        f.write(f"  mol pqr {output_pqr}\n")
        f.write("end\n")

        f.write("elec name solvated\n")
        f.write("  mg-auto\n")
        f.write("  chgm spl2\n")
        f.write("  dime 289 289 289\n")
        f.write("  cglen 150.0 150.0 150.0\n")
        f.write("  fglen 150.0 150.0 150.0\n")
        f.write("  cgcent mol 1\n")
        f.write("  fgcent mol 1\n")
        f.write("  mol 1\n")
        f.write("  lpbe\n")
        f.write("  bcfl mdh\n")
        f.write("  pdie 2.0\n")
        f.write("  sdie 78.54\n")
        f.write("  sdens 10.0\n")
        f.write("  srad 1.4\n")
        f.write("  swin 0.3\n")
        f.write("  srfm smol\n")
        f.write("  temp 298.15\n")
        f.write("  calcenergy total\n")
        f.write("  calcforce no\n")
        f.write(f"  write pot dx {output_pqr}_solv\n")
        f.write("end\n")

        f.write("elec name vacuum\n")
        f.write("  mg-auto\n")
        f.write("  chgm spl2\n")
        f.write("  dime 289 289 289\n")
        f.write("  cglen 150.0 150.0 150.0\n")
        f.write("  fglen 150.0 150.0 150.0\n")
        f.write("  cgcent mol 1\n")
        f.write("  fgcent mol 1\n")
        f.write("  mol 1\n")
        f.write("  lpbe\n")
        f.write("  bcfl sdh\n")
        f.write("  pdie 2.0\n")
        f.write("  sdie 1\n")
        f.write("  sdens 10.0\n")
        f.write("  srad 1.4\n")
        f.write("  swin 0.3\n")
        f.write("  srfm smol\n")
        f.write("  temp 298.15\n")
        f.write("  calcenergy total\n")
        f.write("  calcforce no\n")
        f.write(f"  write pot dx {output_pqr}_vac\n")
        f.write("end\n")

        f.write("print elecEnergy solvated - vacuum end\n")
        f.write("quit\n")

# function for running apbs
def run_apbs(output_apbs_in, apbs_out):
    command = ["apbs", output_apbs_in]
    with open(apbs_out, "w") as output_file:
        subprocess.run(command, stdout=output_file)

# function for extracting energy value
def extract_global_energy(apbs_out):
    with open(apbs_out, 'r') as apbs_output_file:
        lines = apbs_output_file.readlines()
        for line in lines:
            if "Global net ELEC energy" in line:
                global_energy_value_sci = line.split("=")[1].strip().split()[0]
                global_energy_value = "{:.6f}".format(float(global_energy_value_sci))
                return global_energy_value

# function for electrostatics calculations and export
def process_single_pdb(file_path, root_directory, output_file):
    helix_content, sheet_content, chain_length = calculate_secondary_structure(file_path)
    structured_content = helix_content + sheet_content
    relative_path = os.path.relpath(file_path, root_directory)

    protein_name = os.path.basename(os.path.dirname(file_path))
    output_folder = f"/project/target/complete/{protein_name}/apbs/"
    os.makedirs(output_folder, exist_ok=True)

    output_pqr = f"{output_folder}{protein_name}.pqr"
    output_apbs_in = f"{output_folder}{protein_name}_apbs_input.in"
    apbs_out = f"{output_folder}{protein_name}_apbs.out"

    if os.path.exists(apbs_out):
        print(f"\033[0;35m[INFO] Poisson-Boltzmann electrostatics calculation for {protein_name} already exists, skipping.\033[0m")
        solvation_electrostatic_energy = extract_global_energy(apbs_out)
        return (relative_path, helix_content, sheet_content, structured_content, chain_length, output_pqr, output_apbs_in, solvation_electrostatic_energy)

    print(f"\033[0;35m[INFO] Calculating electrostatics for {protein_name}.\033[0m")
    run_pdb2pqr(file_path, output_pqr)
    generate_apbs_input(output_pqr, output_apbs_in)
    run_apbs(output_apbs_in, apbs_out)
    solvation_electrostatic_energy = extract_global_energy(apbs_out)
    return (relative_path, helix_content, sheet_content, structured_content, chain_length, output_pqr, output_apbs_in, solvation_electrostatic_energy)

# function for reading files in parallel
def process_pdb_files_parallel(root_directory, output_file):
    results = []
    max_workers = max(1, os.cpu_count() // 2)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for subdir, dirs, files in os.walk(root_directory):
            for pdb_file in files:
                if pdb_file.endswith(".pdb") and "unrelaxed" not in pdb_file and "ranked" not in pdb_file:
                    file_path = os.path.join(subdir, pdb_file)
                    future = executor.submit(process_single_pdb, file_path, root_directory, output_file)
                    futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)

    return results

# function for writing results
def save_results_to_file(results, output_file):
    with open(output_file, "w") as f:
        f.write("filename\talpha\tbeta\t%structured\tchain_length\tsolvation_electrostatic_energy_kJ/mol\n")
        sorted_results = natsorted(results, key=lambda x: x[0])
        for result in sorted_results:
            f.write(f"{result[0]}\t{result[1]:.2f}\t{result[2]:.2f}\t{result[3]:.2f}\t{result[4]}\t{result[7]}\n")

# execute
if __name__ == "__main__":
    input_directory = "/project/target/complete"
    output_file = "/project/target/ss_summary.tsv"

    processed_results = process_pdb_files_parallel(input_directory, output_file)
    save_results_to_file(processed_results, output_file)
    print(f"\033[0;35m[INFO] Secondary structure percentages, chain lengths, and electrostatics saved to {output_file}\033[0m")