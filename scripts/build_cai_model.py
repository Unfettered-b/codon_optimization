#!/usr/bin/env python3

import json
from collections import defaultdict
from Bio import SeqIO

#############################################
# Inputs / Outputs
#############################################

cds_fasta = snakemake.input.cds
output_json = snakemake.output[0]

# Optional ID file
ids_file = snakemake.input.get("ids", None)

#############################################
# Genetic Code (no stops)
#############################################

aa_to_codons = {
    "F": ["TTT","TTC"],
    "L": ["TTA","TTG","CTT","CTC","CTA","CTG"],
    "I": ["ATT","ATC","ATA"],
    "M": ["ATG"],
    "V": ["GTT","GTC","GTA","GTG"],
    "S": ["TCT","TCC","TCA","TCG","AGT","AGC"],
    "P": ["CCT","CCC","CCA","CCG"],
    "T": ["ACT","ACC","ACA","ACG"],
    "A": ["GCT","GCC","GCA","GCG"],
    "Y": ["TAT","TAC"],
    "H": ["CAT","CAC"],
    "Q": ["CAA","CAG"],
    "N": ["AAT","AAC"],
    "K": ["AAA","AAG"],
    "D": ["GAT","GAC"],
    "E": ["GAA","GAG"],
    "C": ["TGT","TGC"],
    "W": ["TGG"],
    "R": ["CGT","CGC","CGA","CGG","AGA","AGG"],
    "G": ["GGT","GGC","GGA","GGG"]
}

codon_to_aa = {}
for aa, codons in aa_to_codons.items():
    for c in codons:
        codon_to_aa[c] = aa

#############################################
# Load optional ID subset
#############################################

if ids_file:
    with open(ids_file) as f:
        allowed_ids = set(line.strip() for line in f)
else:
    allowed_ids = None

#############################################
# Count codons
#############################################

codon_counts = defaultdict(int)

for record in SeqIO.parse(cds_fasta, "fasta"):

    if allowed_ids is not None and record.id not in allowed_ids:
        continue

    seq = str(record.seq).upper()

    # Ensure full codons only
    seq = seq[:len(seq) - (len(seq) % 3)]

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in codon_to_aa:
            codon_counts[codon] += 1

#############################################
# Compute CAI weights with Laplace smoothing
#############################################

alpha = 1  # pseudocount

codon_weights = {}

for aa, codons in aa_to_codons.items():

    k = len(codons)
    total = sum(codon_counts[c] for c in codons)

    # Apply Laplace smoothing
    smoothed_freqs = {
        c: (codon_counts[c] + alpha) / (total + alpha * k)
        for c in codons
    }

    max_freq = max(smoothed_freqs.values())

    for c in codons:
        codon_weights[c] = smoothed_freqs[c] / max_freq

#############################################
# Write JSON
#############################################

with open(output_json, "w") as out:
    json.dump(codon_weights, out, indent=4)