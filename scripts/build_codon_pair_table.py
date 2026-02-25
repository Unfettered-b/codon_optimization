#!/usr/bin/env python3
import json
import math
import collections
from Bio import SeqIO

cds_fasta = snakemake.input.cds
out_path  = snakemake.output.pairs

pair_counts   = collections.Counter()
single_counts = collections.Counter()

for rec in SeqIO.parse(cds_fasta, "fasta"):
    seq = str(rec.seq).upper().replace("U", "T")
    if len(seq) % 3 != 0:
        continue
    codons = [seq[i: i + 3] for i in range(0, len(seq), 3)]
    for c in codons:
        single_counts[c] += 1
    for a, b in zip(codons, codons[1:]):
        pair_counts[a + b] += 1

total_singles = sum(single_counts.values())
total_pairs   = sum(pair_counts.values())

scores = {}
for pair, obs in pair_counts.items():
    a, b     = pair[:3], pair[3:]
    expected = (
        (single_counts[a] / total_singles)
        * (single_counts[b] / total_singles)
        * total_pairs
    )
    if expected > 0:
        scores[pair] = math.log(obs / expected)

with open(out_path, "w") as f:
    json.dump(scores, f, indent=2)