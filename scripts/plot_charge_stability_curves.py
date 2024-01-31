#!/usr/bin/env python3

# This script reads the output of propka   #
# and creates plots of the stability and   #
# folded and unfolded charge states of     #
# proteins.                                #

# Steven Heaton
# 2024-02-01

import os
import pandas as pd
import matplotlib.pyplot as plt

# extract data from propka output
def extract_charges(charge_file_path):

    df_charge = pd.read_csv(charge_file_path, delimiter='\t')

    pH_values_charge = df_charge['pH']
    unfolded_values = df_charge['unfolded']
    folded_values = df_charge['folded']
    return pH_values_charge, unfolded_values, folded_values

def extract_stability(stability_file_path):

    df_stability = pd.read_csv(stability_file_path, delimiter='\t')

    pH_values_stability = df_stability['pH']
    free_energy = df_stability['free_energy']
    return pH_values_stability, free_energy

# function for plotting charge curves and annotations
def plot_charge_curves(charge_file_path):

    pH_values_charge, unfolded_values, folded_values = extract_charges(charge_file_path)

    min_y = min(pd.concat([unfolded_values, folded_values]))
    max_y = max(pd.concat([unfolded_values, folded_values]))
    plt.ylim(min_y, max_y)

    filename_extensionless = os.path.splitext(os.path.basename(charge_file_path))[0]
    protein_name = filename_extensionless.split('_')[0]
    plt.title(f"{protein_name}: Charge vs. pH")

    # pH/charge values
    plt.plot(pH_values_charge, folded_values, label='folded', color='blue')
    plt.plot(pH_values_charge, unfolded_values, label='unfolded', color='red')
    plt.xlabel('pH')
    plt.ylabel('Charge')
    
    # define areas of interest
    unfolded_pi = pH_values_charge.iloc[(unfolded_values - 0).abs().argsort()[:1]].values[0]
    folded_pi = pH_values_charge.iloc[(folded_values - 0).abs().argsort()[:1]].values[0]
    charge_at7_unfolded = unfolded_values[pH_values_charge.sub(7).abs().idxmin()]
    charge_at7_folded = folded_values[pH_values_charge.sub(7).abs().idxmin()]

    plt.axvline(x=unfolded_pi, linestyle='--', color='red', label=f'pI (unfolded)={unfolded_pi}')
    plt.axvline(x=folded_pi, linestyle='--', color='blue', label=f'pI (folded)={folded_pi}')
    plt.axvline(x=7, linestyle='dotted', color='red', label=f'charge (pH7, unfolded)={charge_at7_unfolded:.2f}')
    plt.axvline(x=7, linestyle='dotted', color='blue', label=f'charge (pH7, folded)={charge_at7_folded:.2f}')

    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')

    charge_output_filename = f"{protein_name}_charge_plot.png"
    charge_output_path = os.path.join(os.path.dirname(charge_file_path), charge_output_filename)
    plt.savefig(charge_output_path, dpi=600, bbox_inches='tight')

    plt.clf()

# function for plotting stability curves and annotations
def plot_stability_curves(stability_file_path):
    pH_values_stability, free_energy = extract_stability(stability_file_path)

    filename_extensionless = os.path.splitext(os.path.basename(stability_file_path))[0]
    protein_name = filename_extensionless.split('_')[0]
    plt.title(f"{protein_name}: Free energy of folding vs. pH")

    # pH/stability values
    plt.plot(pH_values_stability, free_energy, color='green')
    plt.xlabel('pH')
    plt.ylabel('Free energy of folding (kcal/mol)')

    # define areas of interest
    min_idx = free_energy.idxmin()
    min_pH_values = pH_values_stability[free_energy == free_energy.iloc[min_idx]]
    plt.fill_between(min_pH_values, free_energy.min(), free_energy.max(), color='gray', alpha=0.3, label=f'pH range of max stability')

    if min_pH_values.nunique() == 1:
        legend_text = f'pH of max. stability: {min_pH_values.iloc[0]}'
    else:
        legend_text = f'pH range of max stability: {min_pH_values.min()} - {min_pH_values.max()}'

    plt.legend([legend_text], loc='upper left', bbox_to_anchor=(1.05, 0.5))

    stability_output_filename = f"{protein_name}_stability_plot.png"
    stability_output_path = os.path.join(os.path.dirname(stability_file_path), stability_output_filename)
    plt.savefig(stability_output_path, dpi=600, bbox_inches='tight')

    plt.clf()


# os.walk through subdirectories
root_directory = "/project/target/complete"
for foldername, subfolders, filenames in os.walk(root_directory):
    for filename in filenames:
        if filename.endswith("_charge_values.tsv"):
            charge_file_path = os.path.join(foldername, filename)
            plot_charge_curves(charge_file_path)
        if filename.endswith("_stability_values.tsv"):
            stability_file_path = os.path.join(foldername, filename)
            plot_stability_curves(stability_file_path)
