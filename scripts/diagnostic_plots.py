#!/usr/bin/env python3

import json
import math
import os
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

#############################################
# Snakemake I/O
#############################################

all_cds_fasta = snakemake.input.all_cds
gene_metrics_tsv = snakemake.input.gene_metrics
ribosomal_ids_file = snakemake.input.ribosomal_ids
ribo_weights_file = snakemake.input.ribo_weights
genome_weights_file = snakemake.input.genome_weights
final_seq_fasta = snakemake.input.final_seq
pareto_front_tsv = snakemake.input.pareto_front
pareto_sequences_fasta = snakemake.input.pareto_sequences
ga_log_tsv = snakemake.input.ga_log
output_dir = snakemake.output[0]

os.makedirs(output_dir, exist_ok=True)

#############################################
# Constants / helpers
#############################################

AA_TO_CODONS = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
}

SENSE_CODONS = [c for codons in AA_TO_CODONS.values() for c in codons]


def codons_from_seq(seq):
    seq = str(seq).upper()
    n = len(seq) - (len(seq) % 3)
    return [seq[i : i + 3] for i in range(0, n, 3)]


def compute_gc(seq):
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_cai(seq, weights):
    codons = codons_from_seq(seq)
    if not codons:
        return 0.0
    logs = [math.log(max(1e-12, weights.get(c, 1e-8))) for c in codons]
    return math.exp(sum(logs) / len(logs))


def codon_freq_vector(seq):
    codons = codons_from_seq(seq)
    counts = Counter(codons)
    total = sum(counts.values())
    if total == 0:
        return np.zeros(len(SENSE_CODONS))
    return np.array([counts.get(c, 0) / total for c in SENSE_CODONS], dtype=float)


def rscu_from_counts(counts):
    rscu = {}
    for aa, codons in AA_TO_CODONS.items():
        family_total = sum(counts.get(c, 0) for c in codons)
        if family_total == 0:
            for c in codons:
                rscu[c] = 0.0
            continue
        expected = family_total / len(codons)
        for c in codons:
            rscu[c] = counts.get(c, 0) / expected if expected > 0 else 0.0
    return rscu


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=250)
    plt.close()


#############################################
# Load data
#############################################

metrics = pd.read_csv(gene_metrics_tsv, sep="\t")
ga_log = pd.read_csv(ga_log_tsv, sep="\t")
pareto = pd.read_csv(pareto_front_tsv, sep="\t") if os.path.exists(pareto_front_tsv) else pd.DataFrame()

with open(ribo_weights_file) as f:
    ribo_weights = json.load(f)
with open(genome_weights_file) as f:
    genome_weights = json.load(f)

with open(ribosomal_ids_file) as f:
    ribosomal_ids = {line.strip() for line in f if line.strip()}

cds_records = list(SeqIO.parse(all_cds_fasta, "fasta"))
final_record = next(SeqIO.parse(final_seq_fasta, "fasta"))
final_seq = str(final_record.seq).upper()

pareto_records = list(SeqIO.parse(pareto_sequences_fasta, "fasta")) if os.path.exists(pareto_sequences_fasta) else []

# Partition CDS
all_cds_seqs = [str(r.seq).upper() for r in cds_records]
ribo_cds_seqs = [str(r.seq).upper() for r in cds_records if r.id in ribosomal_ids]
if not ribo_cds_seqs:
    ribo_cds_seqs = all_cds_seqs[:]

# Shared targets
gc_target = float(ga_log["gc_target"].iloc[-1]) if "gc_target" in ga_log.columns else compute_gc("".join(all_cds_seqs)[:10000])
sigma = float(snakemake.config.get("gc", {}).get("sigma", 0.05))

#############################################
# 1 & 2 ENC/GC3 histograms
#############################################

plt.figure(figsize=(6, 4))
plt.hist(metrics["ENC"].dropna(), bins=40, alpha=0.8, color="#4472c4", edgecolor="black")
plt.axvline(metrics["ENC"].median(), color="red", linestyle="--", label=f"Median={metrics['ENC'].median():.1f}")
plt.xlabel("ENC")
plt.ylabel("Gene count")
plt.title("ENC distribution across CDS")
plt.legend()
savefig(os.path.join(output_dir, "enc_distribution.png"))

plt.figure(figsize=(6, 4))
plt.hist(metrics["GC3"].dropna(), bins=40, alpha=0.8, color="#70ad47", edgecolor="black")
plt.axvline(metrics["GC3"].median(), color="red", linestyle="--", label=f"Median={metrics['GC3'].median():.2f}")
plt.xlabel("GC3")
plt.ylabel("Gene count")
plt.title("GC3 distribution across CDS")
plt.legend()
savefig(os.path.join(output_dir, "gc3_distribution.png"))

#############################################
# 3 Wright hexbin + residual
#############################################

x = metrics["GC3"].to_numpy()
y = metrics["ENC"].to_numpy()
mask = np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]

plt.figure(figsize=(6, 5))
if len(x) > 0:
    plt.hexbin(x, y, gridsize=35, cmap="viridis", mincnt=1)
    cb = plt.colorbar()
    cb.set_label("Gene density")
curve_x = np.linspace(0.01, 0.99, 200)
curve_y = 2 + curve_x + 29 / (curve_x**2 + (1 - curve_x) ** 2)
plt.plot(curve_x, curve_y, color="red", linewidth=1.5, label="Wright expected")
plt.xlabel("GC3")
plt.ylabel("ENC")
plt.title("Wright density plot")
plt.legend()
savefig(os.path.join(output_dir, "wright_hexbin.png"))

expected = 2 + x + 29 / (x**2 + (1 - x) ** 2)
resid = y - expected
plt.figure(figsize=(6, 4))
plt.axhline(0, color="black", linewidth=1)
plt.scatter(x, resid, s=10, alpha=0.4)
plt.xlabel("GC3")
plt.ylabel("ENC - Wright expected")
plt.title("Wright residuals")
savefig(os.path.join(output_dir, "wright_residuals.png"))

#############################################
# 4 RSCU heatmap (genome / ribosomal / final)
#############################################

def counts_from_seqs(seqs):
    cnt = Counter()
    for s in seqs:
        cnt.update(codons_from_seq(s))
    return cnt

cnt_genome = counts_from_seqs(all_cds_seqs)
cnt_ribo = counts_from_seqs(ribo_cds_seqs)
cnt_final = Counter(codons_from_seq(final_seq))

rscu_genome = rscu_from_counts(cnt_genome)
rscu_ribo = rscu_from_counts(cnt_ribo)
rscu_final = rscu_from_counts(cnt_final)

rows = []
row_labels = []
for aa, codons in AA_TO_CODONS.items():
    for c in codons:
        rows.append([rscu_genome.get(c, 0), rscu_ribo.get(c, 0), rscu_final.get(c, 0)])
        row_labels.append(f"{aa}:{c}")
mat = np.array(rows)

plt.figure(figsize=(7, 12))
plt.imshow(mat, aspect="auto", cmap="magma")
plt.colorbar(label="RSCU")
plt.xticks([0, 1, 2], ["Genome", "Ribosomal", "Final"])
plt.yticks(range(len(row_labels)), row_labels, fontsize=5)
plt.title("RSCU heatmap")
savefig(os.path.join(output_dir, "rscu_heatmap.png"))

#############################################
# 5 CAI distributions with final marker
#############################################

cai_genome_w = [compute_cai(s, genome_weights) for s in all_cds_seqs]
cai_ribo_w = [compute_cai(s, ribo_weights) for s in all_cds_seqs]
final_cai_genome = compute_cai(final_seq, genome_weights)
final_cai_ribo = compute_cai(final_seq, ribo_weights)

plt.figure(figsize=(7, 4))
plt.hist(cai_genome_w, bins=40, alpha=0.6, label="CDS CAI (genome weights)")
plt.hist(cai_ribo_w, bins=40, alpha=0.4, label="CDS CAI (ribosomal weights)")
plt.axvline(final_cai_genome, color="red", linestyle="--", label=f"Final CAI genome={final_cai_genome:.3f}")
plt.axvline(final_cai_ribo, color="purple", linestyle=":", label=f"Final CAI ribo={final_cai_ribo:.3f}")
plt.xlabel("CAI")
plt.ylabel("Gene count")
plt.title("CAI distributions and optimized sequence")
plt.legend(fontsize=8)
savefig(os.path.join(output_dir, "cai_distributions.png"))

#############################################
# 6 Sliding-window GC profile
#############################################

win = int(snakemake.config.get("fungal", {}).get("gc_window", 50))
vals = []
positions = []
for i in range(0, max(1, len(final_seq) - win + 1)):
    w = final_seq[i : i + win]
    if len(w) < win:
        continue
    positions.append(i)
    vals.append(compute_gc(w))

plt.figure(figsize=(9, 3))
if vals:
    plt.plot(positions, vals, linewidth=1)
plt.axhline(gc_target, color="green", linestyle="--", label=f"Target={gc_target:.3f}")
plt.axhline(gc_target + sigma, color="orange", linestyle=":", label=f"±1σ ({sigma:.2f})")
plt.axhline(gc_target - sigma, color="orange", linestyle=":")
plt.ylim(0, 1)
plt.xlabel("Nucleotide position")
plt.ylabel(f"GC ({win} nt window)")
plt.title("Sliding-window GC profile (final sequence)")
plt.legend(fontsize=8)
savefig(os.path.join(output_dir, "sliding_gc_profile.png"))

#############################################
# 7 Rare codon run-length diagnostics
#############################################

final_codons = codons_from_seq(final_seq)
rare = [1 if genome_weights.get(c, 1e-8) < 0.1 else 0 for c in final_codons]
runs = []
cur = 0
for v in rare:
    if v == 1:
        cur += 1
    elif cur > 0:
        runs.append(cur)
        cur = 0
if cur > 0:
    runs.append(cur)

plt.figure(figsize=(6, 4))
if runs:
    bins = np.arange(1, max(runs) + 2)
    plt.hist(runs, bins=bins, align="left", rwidth=0.8)
else:
    plt.text(0.5, 0.5, "No rare-codon runs detected", ha="center", va="center")
plt.xlabel("Rare-codon run length")
plt.ylabel("Count")
plt.title("Rare codon run-length distribution")
savefig(os.path.join(output_dir, "rare_codon_runs.png"))

#############################################
# 8 Pareto scatter matrix (or fallback)
#############################################

pareto_cols = [c for c in ["cai_r", "cai_g", "gc_score", "splice_score"] if c in pareto.columns]
if len(pareto) >= 2 and len(pareto_cols) >= 2:
    pd.plotting.scatter_matrix(pareto[pareto_cols], figsize=(8, 8), alpha=0.7, diagonal="hist")
    plt.suptitle("Pareto objective scatter matrix", y=1.02)
    savefig(os.path.join(output_dir, "pareto_scatter_matrix.png"))
else:
    plt.figure(figsize=(6, 3))
    plt.text(0.5, 0.5, "Pareto scatter unavailable\n(needs >=2 candidates)", ha="center", va="center")
    plt.axis("off")
    savefig(os.path.join(output_dir, "pareto_scatter_matrix.png"))

#############################################
# 9 Candidate diversity diagnostics
#############################################

seqs = [str(r.seq).upper() for r in pareto_records]
pairwise = []
for i in range(len(seqs)):
    for j in range(i + 1, len(seqs)):
        a, b = seqs[i], seqs[j]
        n = min(len(a), len(b))
        if n == 0:
            continue
        dist = sum(1 for k in range(0, n, 3) if a[k : k + 3] != b[k : k + 3])
        pairwise.append(dist / max(1, n // 3))

plt.figure(figsize=(6, 4))
if pairwise:
    plt.hist(pairwise, bins=20)
    plt.xlabel("Pairwise codon Hamming distance (fraction)")
    plt.ylabel("Pair count")
else:
    plt.text(0.5, 0.5, "Insufficient candidates for diversity histogram", ha="center", va="center")
plt.title("Top-candidate diversity")
savefig(os.path.join(output_dir, "candidate_diversity.png"))

#############################################
# 10 Objective trajectory (GA diagnostics)
#############################################

plt.figure(figsize=(8, 4))
plt.plot(ga_log["generation"], ga_log["best_fitness"], label="Best fitness")
plt.plot(ga_log["generation"], ga_log["mean_fitness"], label="Mean fitness")
if "pareto_front_size" in ga_log.columns:
    plt.plot(ga_log["generation"], ga_log["pareto_front_size"], label="Pareto front size")
if "mutation_rate" in ga_log.columns:
    plt.plot(ga_log["generation"], ga_log["mutation_rate"], label="Mutation rate")
plt.xlabel("Generation")
plt.ylabel("Metric value")
plt.title("GA objective / control trajectories")
plt.legend(fontsize=8)
savefig(os.path.join(output_dir, "ga_objective_trajectories.png"))

#############################################
# 11 Synonymous entropy (genome/ribo/final)
#############################################

def synonymous_entropy(counts):
    out = {}
    for aa, codons in AA_TO_CODONS.items():
        vals = np.array([counts.get(c, 0) for c in codons], dtype=float)
        if vals.sum() == 0:
            out[aa] = 0.0
            continue
        p = vals / vals.sum()
        nz = p[p > 0]
        h = -np.sum(nz * np.log2(nz))
        out[aa] = h / np.log2(len(codons)) if len(codons) > 1 else 0.0
    return out

h_gen = synonymous_entropy(cnt_genome)
h_rib = synonymous_entropy(cnt_ribo)
h_fin = synonymous_entropy(cnt_final)

aas = list(AA_TO_CODONS.keys())
xi = np.arange(len(aas))
width = 0.27

plt.figure(figsize=(10, 4))
plt.bar(xi - width, [h_gen[a] for a in aas], width=width, label="Genome")
plt.bar(xi, [h_rib[a] for a in aas], width=width, label="Ribosomal")
plt.bar(xi + width, [h_fin[a] for a in aas], width=width, label="Final")
plt.xticks(xi, aas)
plt.ylim(0, 1.05)
plt.ylabel("Normalized synonymous entropy")
plt.title("Synonymous codon entropy by amino acid")
plt.legend(fontsize=8)
savefig(os.path.join(output_dir, "synonymous_entropy.png"))

#############################################
# 12 PCA of codon usage vectors
#############################################

X = np.array([codon_freq_vector(s) for s in all_cds_seqs])
if len(X) >= 2:
    mu = X.mean(axis=0, keepdims=True)
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pcs = Xc @ Vt.T[:, :2]
    final_pc = (codon_freq_vector(final_seq) - mu.ravel()) @ Vt.T[:, :2]

    plt.figure(figsize=(6, 5))
    plt.scatter(pcs[:, 0], pcs[:, 1], s=8, alpha=0.35, label="Genome CDS")
    plt.scatter(final_pc[0], final_pc[1], color="red", s=60, marker="*", label="Final optimized")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of codon-usage vectors")
    plt.legend(fontsize=8)
    savefig(os.path.join(output_dir, "codon_usage_pca.png"))
else:
    plt.figure(figsize=(6, 3))
    plt.text(0.5, 0.5, "Insufficient CDS vectors for PCA", ha="center", va="center")
    plt.axis("off")
    savefig(os.path.join(output_dir, "codon_usage_pca.png"))

#############################################
# 13 Motif risk map
#############################################

fungal_cfg = snakemake.config.get("fungal", {})
splice_cfg = snakemake.config.get("splice", {})
polya_motifs = fungal_cfg.get("polya_motifs", ["AATAAA", "ATTAAA"])
donor_motifs = fungal_cfg.get("donor_motifs", ["AGGT", "CAGGT", "AAGGT", "GTATGT"])
min_dist = int(splice_cfg.get("intron_min_distance", 20))
max_dist = int(splice_cfg.get("intron_max_distance", 200))


def find_hits(seq, motifs):
    hits = []
    for m in motifs:
        start = 0
        while True:
            pos = seq.find(m, start)
            if pos == -1:
                break
            hits.append((m, pos))
            start = pos + 1
    return hits


def find_gtag_pairs(seq):
    pairs = []
    for i in range(len(seq) - 2):
        if seq[i : i + 2] == "GT":
            window = seq[i + min_dist : i + max_dist]
            rel = window.find("AG")
            if rel != -1:
                pairs.append((i, i + min_dist + rel))
    return pairs

polya_hits = find_hits(final_seq, polya_motifs)
donor_hits = find_hits(final_seq, donor_motifs)
pairs = find_gtag_pairs(final_seq)

plt.figure(figsize=(10, 2.8))
plt.hlines(0.5, 0, len(final_seq), color="lightgray", linewidth=2)
for _, pos in polya_hits:
    plt.vlines(pos, 0.62, 0.88, color="purple", linewidth=1)
for _, pos in donor_hits:
    plt.vlines(pos, 0.12, 0.38, color="red", linewidth=1)
for s, e in pairs:
    plt.plot([s, e], [0.5, 0.5], color="orange", alpha=0.5, linewidth=1)
plt.text(0, 0.9, "PolyA motifs", color="purple", fontsize=8)
plt.text(0, 0.08, "Donor motifs", color="red", fontsize=8)
plt.text(0, 0.53, "GT...AG spans", color="orange", fontsize=8)
plt.xlim(0, max(1, len(final_seq)))
plt.ylim(0, 1)
plt.yticks([])
plt.xlabel("Nucleotide position")
plt.title("Transcript-risk motif map (final sequence)")
savefig(os.path.join(output_dir, "motif_risk_map.png"))
