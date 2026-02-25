#!/usr/bin/env python3

import os
import json
import math
import matplotlib.pyplot as plt
from Bio import SeqIO

#############################################
# Snakemake I/O
#############################################

snapshots_fasta = snakemake.input.snapshots
genome_weights_file = snakemake.input.genome_weights
output_dir = snakemake.output[0]

os.makedirs(output_dir, exist_ok=True)

#############################################
# Load genome weights
#############################################

with open(genome_weights_file) as f:
    genome_weights = json.load(f)

#############################################
# Helper Functions
#############################################

def compute_cai_per_codon(seq, weights):
    """
    Returns list of w values per codon
    """
    w_values = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        w = weights.get(codon, 1e-8)
        w_values.append(w)
    return w_values

#############################################
# Process Each Snapshot
#############################################

for record in SeqIO.parse(snapshots_fasta, "fasta"):

    seq = str(record.seq)
    codon_count = len(seq) // 3

    w_values = compute_cai_per_codon(seq, genome_weights)

    x_positions = list(range(codon_count))

    plt.figure(figsize=(12, 2))

    # Cream backbone block
    plt.axhspan(0, 1, xmin=0, xmax=1, color="#f5f5dc")

    # Plot rare codons
    for i, w in enumerate(w_values):
        if w < 0.1:
            plt.vlines(i, 0, 1, colors="red", linewidth=1)
        elif w < 0.2:
            plt.vlines(i, 0, 1, colors="blue", linewidth=1)

    plt.xlim(0, codon_count)
    plt.ylim(0, 1)

    plt.xticks([])
    plt.yticks([])

    plt.title(f"{record.id} — CAI Codon Usage")

    output_path = os.path.join(output_dir, f"{record.id}.png")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()