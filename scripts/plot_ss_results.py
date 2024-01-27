#!/usr/bin/env python3

# This script plots a 2D scatter plot #
# from the length/structure/charge    #
# calculations in ss_summary.tsv      #

# Steven Heaton
# 2024-01-18

import pandas as pd
import matplotlib.pyplot as plt

# read data and extract relevant information
df = pd.read_csv("/project/target/ss_summary.tsv", sep='\t')

chain_length = df['chain_length']
percent_structured = df['%structured']
surface_charge = df['solvation_electrostatic_energy_kJ/mol']
filenames = df['filename']

protein_names = filenames.apply(lambda x: x.split("/")[0])

# plot
label_fontsize = 14
title_fontsize = 18
dims = (32, 28)
dims = (dims[0] / 2.54, dims[1] / 2.54)

fig, ax = plt.subplots(figsize=dims)
scatter = ax.scatter(chain_length, percent_structured, c=surface_charge, cmap=plt.cm.viridis_r, marker='o', s=60)

for i, txt in enumerate(protein_names):
    ax.annotate(txt, (chain_length.iloc[i], percent_structured.iloc[i]),
                xytext=(15, -15), textcoords='offset points', fontsize=label_fontsize, ha='right')

ax.set_xlabel('Chain Length (aa)', fontsize=label_fontsize)
ax.set_ylabel('% Structured (alpha + beta)', fontsize=label_fontsize)
ax.set_title('Protein Chain Length vs % Structured', fontsize=title_fontsize)

cbar = plt.colorbar(scatter)
cbar.set_label('Solvation energy (kJ/mol)', fontsize=label_fontsize)

plt.tight_layout()

# save plot
plt.savefig("/project/target/scatterplot_sizeVssVcharge.pdf", dpi=600, bbox_inches='tight')
print(f"\033[0;35m[INFO] Chain length vs. secondary structure vs. charge scatterplot saved to /project/target/scatterplot_sizeVssVcharge.pdf\033[0m")