#!/usr/bin/env python3

from Bio import SeqIO

input_fasta = snakemake.input[0]
output_fasta = snakemake.output[0]
output_warnings = snakemake.output[1]

#############################################
# Parameters
#############################################

fungal_cfg = snakemake.config.get("fungal", {})
splice_cfg = snakemake.config.get("splice", {})

WINDOW = int(fungal_cfg.get("gc_window", 50))
LOW_GC = float(fungal_cfg.get("low_gc", 0.30))
HIGH_GC = float(fungal_cfg.get("high_gc", 0.70))

POLYA_MOTIFS = fungal_cfg.get("polya_motifs", ["AATAAA", "ATTAAA"])
DONOR_MOTIFS = fungal_cfg.get("donor_motifs", ["AGGT", "CAGGT", "AAGGT", "GTATGT"])

INTRON_MIN_DISTANCE = int(splice_cfg.get("intron_min_distance", 20))
INTRON_MAX_DISTANCE = int(splice_cfg.get("intron_max_distance", 200))

#############################################
# Utility functions
#############################################

def compute_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def sliding_gc(seq):
    flags = []
    for i in range(len(seq) - WINDOW + 1):
        window = seq[i:i+WINDOW]
        gc = compute_gc(window)
        if gc < LOW_GC or gc > HIGH_GC:
            flags.append((i, round(gc, 3)))
    return flags

def find_motifs(seq, motifs):
    hits = []
    for motif in motifs:
        start = 0
        while True:
            pos = seq.find(motif, start)
            if pos == -1:
                break
            hits.append((motif, pos))
            start = pos + 1
    return hits

def find_gt_ag_pairs(seq, min_dist=20, max_dist=200):
    hits = []
    length = len(seq)
    for i in range(length - 2):
        if seq[i:i+2] == "GT":
            window = seq[i+min_dist:i+max_dist]
            rel = window.find("AG")
            if rel != -1:
                hits.append((i, i+min_dist+rel))
    return hits

#############################################
# Main
#############################################

record = next(SeqIO.parse(input_fasta, "fasta"))
seq = str(record.seq).upper()

warnings = []

# ORF checks
if len(seq) % 3 != 0:
    warnings.append("Sequence length not divisible by 3")

for i in range(0, len(seq), 3):
    codon = seq[i:i+3]
    if codon in ["TAA", "TAG", "TGA"] and i < len(seq) - 3:
        warnings.append(f"Internal stop codon at position {i}")

# PolyA motifs
polya_hits = find_motifs(seq, POLYA_MOTIFS)
for motif, pos in polya_hits:
    warnings.append(f"PolyA motif {motif} at position {pos}")

# Donor motifs
donor_hits = find_motifs(seq, DONOR_MOTIFS)
for motif, pos in donor_hits:
    warnings.append(f"Strong donor motif {motif} at position {pos}")

# Intron-like GT...AG
gt_ag_hits = find_gt_ag_pairs(seq, min_dist=INTRON_MIN_DISTANCE, max_dist=INTRON_MAX_DISTANCE)
for start, end in gt_ag_hits:
    warnings.append(f"Potential intron GT...AG from {start} to {end}")

# Sliding GC
gc_flags = sliding_gc(seq)
for pos, gc in gc_flags:
    warnings.append(f"Extreme GC window at {pos} (GC={gc})")

#############################################
# Write outputs
#############################################

with open(output_fasta, "w") as f:
    f.write(f">{record.id}_final\n{seq}\n")

with open(output_warnings, "w") as w:
    if warnings:
        for line in warnings:
            w.write(line + "\n")
    else:
        w.write("No major fungal expression risks detected.\n")