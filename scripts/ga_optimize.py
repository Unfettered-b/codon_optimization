#!/usr/bin/env python3

import json
import math
import random
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

#############################################
# Snakemake I/O
#############################################

protein_fasta = snakemake.input.protein
ribo_weights_file = snakemake.input.ribo_weights
genome_weights_file = snakemake.input.genome_weights

output_fasta = snakemake.output[0]
output_log = snakemake.output[1]
output_plot = snakemake.output[2]
snapshots_fasta = snakemake.output[3]

config = snakemake.config


# clear snapshots output file

open(snapshots_fasta, "w").close()


#############################################
# Load parameters from config
#############################################

POP_SIZE = config["ga"]["population_size"]
GENERATIONS = config["ga"]["generations"]
MUTATION_RATE = config["ga"]["mutation_rate"]
ELITE_SIZE = config["ga"]["elite_size"]
random.seed(config["ga"]["random_seed"])

alpha = config["fitness"]["alpha"]
beta = config["fitness"]["beta"]

GC_TARGET = config["gc"]["target"]
GC_SIGMA = config["gc"]["sigma"]

DONOR_PENALTY = config["splice"]["donor_penalty"]
INTRON_PENALTY = config["splice"]["intron_penalty"]
INTRON_MIN = config["splice"]["intron_min_distance"]
INTRON_MAX = config["splice"]["intron_max_distance"]


SNAPSHOT_COUNT = config['Snapshot']['count']
SN = max(1, GENERATIONS // SNAPSHOT_COUNT) # getting intervals of snapshots

#############################################
# Load weight models
#############################################

with open(ribo_weights_file) as f:
    ribo_weights = json.load(f)

with open(genome_weights_file) as f:
    genome_weights = json.load(f)

#############################################
# Genetic Code
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

#############################################
# Utility Functions
#############################################

def backtranslate_max(protein_seq):
    return [
        max(aa_to_codons[aa], key=lambda c: ribo_weights.get(c, 1e-8))
        for aa in protein_seq
    ]

def codon_list_to_seq(codon_list):
    return "".join(codon_list)

def compute_cai(dna_seq, weights):
    score = 0.0
    L = 0
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        w = weights.get(codon, 1e-8)
        score += math.log(w)
        L += 1
    return math.exp(score / L)

def compute_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def gc_penalty(seq):
    gc = compute_gc(seq)
    return math.exp(-((gc - GC_TARGET) ** 2) / (2 * GC_SIGMA ** 2))

def splice_penalty(seq):
    penalty = 1.0
    length = len(seq)

    high_risk = ["AGGT", "CAGGT", "AAGGT", "GTATGT"]

    for motif in high_risk:
        if motif in seq:
            penalty *= DONOR_PENALTY

    for i in range(length - 2):
        if seq[i:i+2] == "GT":
            window = seq[i+INTRON_MIN:i+INTRON_MAX]
            if "AG" in window:
                penalty *= INTRON_PENALTY

    return penalty

def fitness(codon_list):
    dna = codon_list_to_seq(codon_list)

    cai_r = compute_cai(dna, ribo_weights)
    cai_g = compute_cai(dna, genome_weights)

    return (
        (cai_r ** alpha)
        * (cai_g ** beta)
        * gc_penalty(dna)
        * splice_penalty(dna)
    )

def mutate(codon_list, protein_seq):
    new = codon_list.copy()
    for i, aa in enumerate(protein_seq):
        if random.random() < MUTATION_RATE:
            new[i] = random.choice(aa_to_codons[aa])
    return new

def tournament(pop, k=3):
    selected = random.sample(pop, k)
    return max(selected, key=lambda ind: ind["fitness"])

#############################################
# Main GA
#############################################

record = next(SeqIO.parse(protein_fasta, "fasta"))
protein_seq = str(record.seq)

base = backtranslate_max(protein_seq)

population = []

for _ in range(POP_SIZE):
    codons = mutate(base, protein_seq)
    population.append({
        "codons": codons,
        "fitness": fitness(codons)
    })

log_rows = []

for gen in range(GENERATIONS):

    population.sort(key=lambda x: x["fitness"], reverse=True)

    best_fit = population[0]["fitness"]
    mean_fit = sum(ind["fitness"] for ind in population) / POP_SIZE

    log_rows.append({
        "generation": gen,
        "best_fitness": best_fit,
        "mean_fitness": mean_fit
    })

    if gen % SN == 0 or gen == GENERATIONS - 1:
        best_seq = codon_list_to_seq(population[0]["codons"])

        with open(snapshots_fasta, "a") as f:
            f.write(f">generation_{gen}_fitness_{best_fit:.6f}\n")
            f.write(best_seq + "\n")
    
    new_population = population[:ELITE_SIZE]

    while len(new_population) < POP_SIZE:
        parent = tournament(population)
        child_codons = mutate(parent["codons"], protein_seq)
        new_population.append({
            "codons": child_codons,
            "fitness": fitness(child_codons)
        })

    population = new_population

#############################################
# Final Best Individual
#############################################

population.sort(key=lambda x: x["fitness"], reverse=True)
best = population[0]
best_seq = codon_list_to_seq(best["codons"])

with open(output_fasta, "w") as f:
    f.write(f">{record.id}_optimized\n{best_seq}\n")

#############################################
# Write GA Log
#############################################

ga_log_df = pd.DataFrame(log_rows)
ga_log_df.to_csv(output_log, sep="\t", index=False)

#############################################
# Fitness vs Generation Plot
#############################################

plt.figure()
plt.plot(ga_log_df["generation"], ga_log_df["best_fitness"], label="Best")
plt.plot(ga_log_df["generation"], ga_log_df["mean_fitness"], label="Mean")
plt.xlabel("Generation")
plt.ylabel("Fitness")
plt.title("Fitness Across Generations")
plt.legend()
plt.tight_layout()
plt.savefig(output_plot)
plt.close()