#!/usr/bin/env python3

# This script reads the output of propka   #
# and creates plots of the folded and      #
# unfolded charge states of proteins.      #

# Steven Heaton
# 2024-01-27

import os
import pandas as pd
import matplotlib.pyplot as plt

# extract data from propka output
def extract_data_from_log(file_path):

    df = pd.read_csv(file_path, delimiter='\t')

    pH_values = df['pH']
    unfolded_values = df['unfolded']
    folded_values = df['folded']
    return pH_values, unfolded_values, folded_values

# function for plotting charge curves and annotations
def plot_charge_curves(file_path):

    pH_values, unfolded_values, folded_values = extract_data_from_log(file_path)

    min_y = min(pd.concat([unfolded_values, folded_values]))
    max_y = max(pd.concat([unfolded_values, folded_values]))
    plt.ylim(min_y, max_y)

    foldername = os.path.basename(os.path.dirname(file_path))
    filename_without_extension = os.path.splitext(os.path.basename(file_path))[0]
    protein_name = filename_without_extension.split('_')[0]
    plt.title(f"{protein_name}: Charge vs. pH")

    # pH/charge values
    plt.plot(pH_values, folded_values, label='folded', color='blue')
    plt.plot(pH_values, unfolded_values, label='unfolded', color='red')
    plt.xlabel('pH')
    plt.ylabel('Charge')
    
    # define areas of interest
    unfolded_pi = pH_values.iloc[(unfolded_values - 0).abs().argsort()[:1]].values[0]
    folded_pi = pH_values.iloc[(folded_values - 0).abs().argsort()[:1]].values[0]
    charge_at7_unfolded = unfolded_values[pH_values.sub(7).abs().idxmin()]
    charge_at7_folded = folded_values[pH_values.sub(7).abs().idxmin()]

    plt.axvline(x=unfolded_pi, linestyle='--', color='red', label=f'pI (unfolded)={unfolded_pi}')
    plt.axvline(x=folded_pi, linestyle='--', color='blue', label=f'pI (folded)={folded_pi}')
    plt.axvline(x=7, linestyle='dotted', color='red', label=f'charge (pH7, unfolded)={charge_at7_unfolded:.2f}')
    plt.axvline(x=7, linestyle='dotted', color='blue', label=f'charge (pH7, folded)={charge_at7_folded:.2f}')

    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')

    output_filename = f"{protein_name}_charge_plot.png"
    output_path = os.path.join(os.path.dirname(file_path), output_filename)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')

    plt.clf()

# os.walk through subdirectories
root_directory = "/project/target/complete"
for foldername, subfolders, filenames in os.walk(root_directory):
    for filename in filenames:
        if filename.endswith(".tsv"):
            file_path = os.path.join(foldername, filename)
            plot_charge_curves(file_path)
