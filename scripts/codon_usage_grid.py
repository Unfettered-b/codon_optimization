#!/usr/bin/env python3

import json
import os
from collections import Counter
import matplotlib.pyplot as plt
from Bio import SeqIO
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import random
import numpy as np

#############################################
# I/O
#############################################

genome_weights_file = snakemake.input.genome_weights
ribo_weights_file = snakemake.input.ribo_weights
final_seq_file = snakemake.input.final_seq
output_dir = snakemake.output[0]

os.makedirs(output_dir, exist_ok=True)

#############################################
# Amino Acid Table
#############################################

aa_table = {
    "F": ("Phenylalanine", "Phe", ["TTT","TTC"]),
    "L": ("Leucine", "Leu", ["TTA","TTG","CTT","CTC","CTA","CTG"]),
    "I": ("Isoleucine", "Ile", ["ATT","ATC","ATA"]),
    "M": ("Methionine", "Met", ["ATG"]),
    "V": ("Valine", "Val", ["GTT","GTC","GTA","GTG"]),
    "S": ("Serine", "Ser", ["TCT","TCC","TCA","TCG","AGT","AGC"]),
    "P": ("Proline", "Pro", ["CCT","CCC","CCA","CCG"]),
    "T": ("Threonine", "Thr", ["ACT","ACC","ACA","ACG"]),
    "A": ("Alanine", "Ala", ["GCT","GCC","GCA","GCG"]),
    "Y": ("Tyrosine", "Tyr", ["TAT","TAC"]),
    "H": ("Histidine", "His", ["CAT","CAC"]),
    "Q": ("Glutamine", "Gln", ["CAA","CAG"]),
    "N": ("Asparagine", "Asn", ["AAT","AAC"]),
    "K": ("Lysine", "Lys", ["AAA","AAG"]),
    "D": ("Aspartate", "Asp", ["GAT","GAC"]),
    "E": ("Glutamate", "Glu", ["GAA","GAG"]),
    "C": ("Cysteine", "Cys", ["TGT","TGC"]),
    "W": ("Tryptophan", "Trp", ["TGG"]),
    "R": ("Arginine", "Arg", ["CGT","CGC","CGA","CGG","AGA","AGG"]),
    "G": ("Glycine", "Gly", ["GGT","GGC","GGA","GGG"])
}

#############################################
# Load Data
#############################################

with open(genome_weights_file) as f:
    genome = json.load(f)

with open(ribo_weights_file) as f:
    ribo = json.load(f)

record = next(SeqIO.parse(final_seq_file, "fasta"))
seq = str(record.seq)
codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
final_counts = Counter(codons)
total = sum(final_counts.values())
final = {k: v/total for k, v in final_counts.items()}

#############################################
# Color Map
#############################################

all_codons = sorted({c for aa in aa_table.values() for c in aa[2]})
palette = list(mcolors.TABLEAU_COLORS.values())
random.seed(42)

codon_colors = {
    codon: palette[i % len(palette)]
    for i, codon in enumerate(all_codons)
}

#############################################
# Pie with Spokes
#############################################

def pie_with_labels(ax, source_dict, codon_list):

    values = np.array([source_dict.get(c, 0) for c in codon_list])
    total = values.sum()

    if total == 0:
        ax.axis("off")
        return

    values = values / total

    wedges, _ = ax.pie(
        values,
        colors=[codon_colors[c] for c in codon_list],
        startangle=90
    )

    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.deg2rad(angle))
        y = np.sin(np.deg2rad(angle))

        ax.annotate(
            f"{codon_list[i]} ({values[i]*100:.1f}%)",
            xy=(x, y),
            xytext=(1.2*x, 1.2*y),
            arrowprops=dict(arrowstyle="-"),
            ha="center",
            va="center",
            fontsize=6
        )

    ax.set_aspect("equal")

#############################################
# Generate Grid (5 residues per image)
#############################################

aa_items = list(aa_table.items())
batch_size = 5

for batch_idx in range(0, len(aa_items), batch_size):

    batch = aa_items[batch_idx:batch_idx+batch_size]
    num_rows = len(batch)

    fig = plt.figure(figsize=(12, num_rows * 2))
    gs = gridspec.GridSpec(num_rows + 1, 4, width_ratios=[1.4,1,1,1])

    # Header row
    headers = ["Residue", "Genome", "Ribosomal", "Final"]
    for col in range(4):
        ax = fig.add_subplot(gs[0, col])
        ax.axis("off")
        ax.text(0.5, 0.5, headers[col], ha="center", va="center", fontsize=11, fontweight="bold")

    # Residue rows
    for row_idx, (aa, (fullname, three, codon_list)) in enumerate(batch, start=1):

        ax_label = fig.add_subplot(gs[row_idx, 0])
        ax_label.axis("off")

        ax_label.text(0.01, 0.7, f"{fullname}", fontsize=10, fontweight="bold")
        ax_label.text(0.01, 0.5, f"{three} ({aa})", fontsize=8)

        for i, codon in enumerate(codon_list):
            ax_label.text(
                0.01 + i*0.07,
                0.3,
                codon,
                color=codon_colors[codon],
                fontsize=8
            )

        pie_with_labels(fig.add_subplot(gs[row_idx,1]), genome, codon_list)
        pie_with_labels(fig.add_subplot(gs[row_idx,2]), ribo, codon_list)
        pie_with_labels(fig.add_subplot(gs[row_idx,3]), final, codon_list)

    plt.tight_layout(pad=1.0)
    plt.savefig(os.path.join(output_dir, f"grid_part_{batch_idx//batch_size + 1}.png"), dpi=300)
    plt.close()