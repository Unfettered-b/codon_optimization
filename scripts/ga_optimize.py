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
genome_cds_fasta = snakemake.input.genome_cds
output_fasta = snakemake.output.optimized
output_log = snakemake.output.log
output_plot = snakemake.output.plot
snapshots_fasta = snakemake.output.snapshots
pareto_front_tsv = snakemake.output.pareto_front
pareto_sequences_fasta = snakemake.output.pareto_sequences
config = snakemake.config

# clear snapshots output file
open(snapshots_fasta, "w").close()

#############################################
# Load parameters from config
#############################################
ga_cfg = config["ga"]
MODE = ga_cfg.get("mode", "scalar").lower()
if MODE not in {"scalar", "pareto"}:
    raise ValueError("ga.mode must be 'scalar' or 'pareto'")

POP_SIZE = int(ga_cfg["population_size"])
GENERATIONS = int(ga_cfg["generations"])
MUTATION_RATE = float(ga_cfg["mutation_rate"])
ELITE_SIZE = int(ga_cfg["elite_size"])
RANDOM_SEED = int(ga_cfg["random_seed"])
PARETO_TOP_K = max(1, int(ga_cfg.get("pareto_top_k", 10)))

# Optional advanced GA controls
CROSSOVER_RATE = float(ga_cfg.get("crossover_rate", 0.75))
TOURNAMENT_SIZE = int(ga_cfg.get("tournament_size", 3))
IMMIGRANT_RATE = float(ga_cfg.get("immigrant_rate", 0.05))
EARLY_STOP_PATIENCE = int(ga_cfg.get("early_stop_patience", 1000))
MIN_IMPROVEMENT = float(ga_cfg.get("min_improvement", 1e-10))
MUTATION_RATE_FINAL = float(ga_cfg.get("mutation_rate_final", MUTATION_RATE))

random.seed(RANDOM_SEED)

alpha = float(config["fitness"]["alpha"])
beta = float(config["fitness"]["beta"])

gc_cfg = config.get("gc", {})
USE_GENOME_GC_TARGET = bool(gc_cfg.get("use_genome_target", True))
GC_TARGET_FALLBACK = float(gc_cfg.get("target", 0.5))
GC_SIGMA = float(gc_cfg.get("sigma", 0.05))

DONOR_PENALTY = float(config["splice"]["donor_penalty"])
INTRON_PENALTY = float(config["splice"]["intron_penalty"])
INTRON_MIN = int(config["splice"]["intron_min_distance"])
INTRON_MAX = int(config["splice"]["intron_max_distance"])

snapshot_cfg = config.get("snapshot", config.get("Snapshot", {}))
SNAPSHOT_COUNT = max(1, int(snapshot_cfg.get("count", 10)))
SN = max(1, GENERATIONS // SNAPSHOT_COUNT)  # snapshot interval

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
    "*": ["TAA", "TAG", "TGA"],  # stop codons — included so the dict never raises KeyError
}

# Reverse lookup: codon -> amino acid (built once at startup for validation)
codon_to_aa = {
    codon: aa
    for aa, codons in aa_to_codons.items()
    for codon in codons
}

#############################################
# Utility Functions
#############################################
def compute_gc(seq):
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_genome_gc_target(cds_fasta_path, fallback):
    total_gc = 0
    total_len = 0
    for rec in SeqIO.parse(cds_fasta_path, "fasta"):
        s = str(rec.seq).upper()
        total_gc += s.count("G") + s.count("C")
        total_len += len(s)
    if total_len == 0:
        return fallback
    return total_gc / total_len


GC_TARGET = (
    compute_genome_gc_target(genome_cds_fasta, GC_TARGET_FALLBACK)
    if USE_GENOME_GC_TARGET
    else GC_TARGET_FALLBACK
)


def backtranslate_max(protein_seq):
    codons = []
    skipped = []
    for pos, aa in enumerate(protein_seq):
        if aa not in aa_to_codons:
            skipped.append((pos, aa))
            continue
        codons.append(max(aa_to_codons[aa], key=lambda c: ribo_weights.get(c, 1e-8)))
    if skipped:
        import warnings
        warnings.warn(
            f"backtranslate_max: skipped {len(skipped)} unrecognised residue(s): "
            + ", ".join(f"pos {p} '{a}'" for p, a in skipped)
        )
    return codons


def codon_list_to_seq(codon_list):
    return "".join(codon_list)


def compute_cai(dna_seq, weights):
    score = 0.0
    L = 0
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i : i + 3]
        if len(codon) != 3:
            continue
        w = weights.get(codon, 1e-8)
        score += math.log(w)
        L += 1
    if L == 0:
        return 0.0
    return math.exp(score / L)


def gc_penalty(seq):
    gc = compute_gc(seq)
    return math.exp(-((gc - GC_TARGET) ** 2) / (2 * GC_SIGMA**2))


def splice_penalty(seq):
    penalty = 1.0
    length = len(seq)
    high_risk = ["AGGT", "CAGGT", "AAGGT", "GTATGT"]
    for motif in high_risk:
        if motif in seq:
            penalty *= DONOR_PENALTY
    for i in range(length - 2):
        if seq[i : i + 2] == "GT":
            window = seq[i + INTRON_MIN : i + INTRON_MAX]
            if "AG" in window:
                penalty *= INTRON_PENALTY
    return penalty


score_cache = {}


def evaluate(codon_list):
    dna = codon_list_to_seq(codon_list)
    cached = score_cache.get(dna)
    if cached is not None:
        return cached
    cai_r = compute_cai(dna, ribo_weights)
    cai_g = compute_cai(dna, genome_weights)
    gc_score = gc_penalty(dna)
    splice_score = splice_penalty(dna)
    scalar_fitness = (cai_r**alpha) * (cai_g**beta) * gc_score * splice_score
    result = {
        "dna": dna,
        "cai_r": cai_r,
        "cai_g": cai_g,
        "gc": compute_gc(dna),
        "gc_score": gc_score,
        "splice_score": splice_score,
        "fitness": scalar_fitness,
    }
    score_cache[dna] = result
    return result


def verify_synonymous(codon_list, protein_seq, context=""):
    """
    Verify that codon_list encodes exactly protein_seq.
    Raises ValueError on mismatch with a detailed diagnostic message.
    Silent and fast when everything is correct.
    """
    if len(codon_list) != len(protein_seq):
        raise ValueError(
            f"[{context}] Length mismatch: {len(codon_list)} codons vs "
            f"{len(protein_seq)} residues."
        )
    mismatches = []
    for i, (codon, aa) in enumerate(zip(codon_list, protein_seq)):
        expected = codon_to_aa.get(codon)
        if expected != aa:
            mismatches.append(
                f"  pos {i}: codon '{codon}' -> '{expected}', expected '{aa}'"
            )
    if mismatches:
        raise ValueError(
            f"[{context}] Non-synonymous mutation(s) detected "
            f"({len(mismatches)} site(s)):\n" + "\n".join(mismatches[:10])
            + ("\n  ..." if len(mismatches) > 10 else "")
        )


def mutate(codon_list, protein_seq, mutation_rate):
    """
    Perform per-codon synonymous mutation.

    Residues absent from aa_to_codons (e.g. ambiguous characters) are silently
    preserved so the codon list length never changes and downstream length
    checks remain valid.
    """
    new = codon_list.copy()
    for i, aa in enumerate(protein_seq):
        # Skip residues not in the codon table — preserve whatever codon is there
        if aa not in aa_to_codons:
            continue
        if random.random() < mutation_rate:
            new[i] = random.choice(aa_to_codons[aa])
    return new


def update_scores(ind, protein_seq=None, debug=False):
    """
    Evaluate an individual and store fitness metrics in-place.

    Parameters
    ----------
    ind : dict
        Individual with a 'codons' key.
    protein_seq : str, optional
        If provided (recommended), run a synonymy check before scoring.
        Pass the global protein_seq so bugs are caught immediately.
    debug : bool
        Set True during development to enable the synonymy check;
        set False (default) for production runs to skip the overhead.
    """
    if debug and protein_seq is not None:
        verify_synonymous(ind["codons"], protein_seq, context="update_scores")
    ind.update(evaluate(ind["codons"]))


#############################################
# GA Operators
#############################################
def crossover(parent_a, parent_b):
    if len(parent_a) <= 1:
        return parent_a.copy(), parent_b.copy()
    cut = random.randint(1, len(parent_a) - 1)
    child_a = parent_a[:cut] + parent_b[cut:]
    child_b = parent_b[:cut] + parent_a[cut:]
    return child_a, child_b


def scheduled_mutation_rate(gen):
    if GENERATIONS <= 1:
        return MUTATION_RATE
    progress = gen / (GENERATIONS - 1)
    return MUTATION_RATE + (MUTATION_RATE_FINAL - MUTATION_RATE) * progress


def dominates(ind_a, ind_b):
    obj_a = (ind_a["cai_r"], ind_a["cai_g"], ind_a["gc_score"], ind_a["splice_score"])
    obj_b = (ind_b["cai_r"], ind_b["cai_g"], ind_b["gc_score"], ind_b["splice_score"])
    not_worse = all(a >= b for a, b in zip(obj_a, obj_b))
    strictly_better = any(a > b for a, b in zip(obj_a, obj_b))
    return not_worse and strictly_better


def fast_non_dominated_sort(pop):
    for p in pop:
        p["dominated_set"] = []
        p["dom_count"] = 0
    fronts = [[]]
    for p in pop:
        for q in pop:
            if p is q:
                continue
            if dominates(p, q):
                p["dominated_set"].append(q)
            elif dominates(q, p):
                p["dom_count"] += 1
        if p["dom_count"] == 0:
            p["rank"] = 0
            fronts[0].append(p)
    i = 0
    while i < len(fronts) and fronts[i]:
        next_front = []
        for p in fronts[i]:
            for q in p["dominated_set"]:
                q["dom_count"] -= 1
                if q["dom_count"] == 0:
                    q["rank"] = i + 1
                    next_front.append(q)
        i += 1
        if next_front:
            fronts.append(next_front)
    return fronts


def assign_crowding_distance(front):
    if not front:
        return
    for ind in front:
        ind["crowding"] = 0.0
    objectives = ["cai_r", "cai_g", "gc_score", "splice_score"]
    n = len(front)
    if n <= 2:
        for ind in front:
            ind["crowding"] = float("inf")
        return
    for obj in objectives:
        front.sort(key=lambda x: x[obj])
        front[0]["crowding"] = float("inf")
        front[-1]["crowding"] = float("inf")
        obj_min = front[0][obj]
        obj_max = front[-1][obj]
        if obj_max == obj_min:
            continue
        for i in range(1, n - 1):
            if math.isinf(front[i]["crowding"]):
                continue
            prev_v = front[i - 1][obj]
            next_v = front[i + 1][obj]
            front[i]["crowding"] += (next_v - prev_v) / (obj_max - obj_min)


def crowded_tournament(pop, k=2):
    contestants = random.sample(pop, min(len(pop), max(2, k)))
    contestants.sort(key=lambda x: (x.get("rank", 10**9), -x.get("crowding", 0.0)))
    return contestants[0]


def scalar_tournament(pop, k=3):
    k = max(2, min(k, len(pop)))
    selected = random.sample(pop, k)
    return max(selected, key=lambda ind: ind["fitness"])


def rank_population(pop):
    fronts = fast_non_dominated_sort(pop)
    for front in fronts:
        assign_crowding_distance(front)
    return fronts


def select_next_generation_nsga2(combined_population, target_size):
    fronts = rank_population(combined_population)
    new_pop = []
    for front in fronts:
        if len(new_pop) + len(front) <= target_size:
            new_pop.extend(front)
        else:
            front.sort(key=lambda x: x.get("crowding", 0.0), reverse=True)
            new_pop.extend(front[: target_size - len(new_pop)])
            break
    return new_pop


#############################################
# Main GA
#############################################
record = next(SeqIO.parse(protein_fasta, "fasta"))
protein_seq = str(record.seq)
base = backtranslate_max(protein_seq)

# Sanity-check the base sequence before the GA starts
verify_synonymous(base, protein_seq, context="backtranslate_max")

population = []
for _ in range(POP_SIZE):
    codons = mutate(base, protein_seq, mutation_rate=MUTATION_RATE)
    ind = {"codons": codons}
    update_scores(ind, protein_seq=protein_seq)
    population.append(ind)

if MODE == "pareto":
    rank_population(population)

log_rows = []
best_overall_fitness = -1.0
stagnation_counter = 0
last_generation = GENERATIONS - 1

for gen in range(GENERATIONS):
    if MODE == "scalar":
        population.sort(key=lambda x: x["fitness"], reverse=True)
    else:
        rank_population(population)
        population.sort(key=lambda x: (x["rank"], -x["crowding"], -x["fitness"]))

    best_fit = population[0]["fitness"]
    mean_fit = sum(ind["fitness"] for ind in population) / POP_SIZE
    front_size = 1
    if MODE == "pareto":
        front_size = sum(1 for ind in population if ind.get("rank", 9999) == 0)

    log_rows.append(
        {
            "generation": gen,
            "best_fitness": best_fit,
            "mean_fitness": mean_fit,
            "mutation_rate": scheduled_mutation_rate(gen),
            "cache_size": len(score_cache),
            "mode": MODE,
            "pareto_front_size": front_size,
            "gc_target": GC_TARGET,
        }
    )

    if best_fit > best_overall_fitness + MIN_IMPROVEMENT:
        best_overall_fitness = best_fit
        stagnation_counter = 0
    else:
        stagnation_counter += 1

    if gen % SN == 0 or gen == GENERATIONS - 1:
        best_seq = population[0]["dna"]
        with open(snapshots_fasta, "a") as f:
            f.write(f">generation_{gen}_fitness_{best_fit:.6f}\n")
            f.write(best_seq + "\n")

    if stagnation_counter >= EARLY_STOP_PATIENCE:
        last_generation = gen
        break

    mutation_rate = scheduled_mutation_rate(gen)

    if MODE == "scalar":
        new_population = population[:ELITE_SIZE]
        while len(new_population) < POP_SIZE:
            parent1 = scalar_tournament(population, k=TOURNAMENT_SIZE)
            parent2 = scalar_tournament(population, k=TOURNAMENT_SIZE)
            if random.random() < CROSSOVER_RATE:
                child1_codons, child2_codons = crossover(parent1["codons"], parent2["codons"])
            else:
                child1_codons = parent1["codons"].copy()
                child2_codons = parent2["codons"].copy()

            child1_codons = mutate(child1_codons, protein_seq, mutation_rate=mutation_rate)
            child1 = {"codons": child1_codons}
            update_scores(child1, protein_seq=protein_seq)
            new_population.append(child1)

            if len(new_population) < POP_SIZE:
                child2_codons = mutate(child2_codons, protein_seq, mutation_rate=mutation_rate)
                child2 = {"codons": child2_codons}
                update_scores(child2, protein_seq=protein_seq)
                new_population.append(child2)

        immigrant_count = max(0, int(POP_SIZE * IMMIGRANT_RATE))
        for _ in range(immigrant_count):
            if len(new_population) <= ELITE_SIZE:
                break
            idx = random.randint(ELITE_SIZE, len(new_population) - 1)
            immigrant = {"codons": mutate(base, protein_seq, mutation_rate=1.0)}
            update_scores(immigrant, protein_seq=protein_seq)
            new_population[idx] = immigrant

        population = new_population

    else:
        rank_population(population)
        offspring = []
        while len(offspring) < POP_SIZE:
            parent1 = crowded_tournament(population, k=TOURNAMENT_SIZE)
            parent2 = crowded_tournament(population, k=TOURNAMENT_SIZE)
            if random.random() < CROSSOVER_RATE:
                child1_codons, child2_codons = crossover(parent1["codons"], parent2["codons"])
            else:
                child1_codons = parent1["codons"].copy()
                child2_codons = parent2["codons"].copy()

            child1_codons = mutate(child1_codons, protein_seq, mutation_rate=mutation_rate)
            child1 = {"codons": child1_codons}
            update_scores(child1, protein_seq=protein_seq)
            offspring.append(child1)

            if len(offspring) < POP_SIZE:
                child2_codons = mutate(child2_codons, protein_seq, mutation_rate=mutation_rate)
                child2 = {"codons": child2_codons}
                update_scores(child2, protein_seq=protein_seq)
                offspring.append(child2)

        immigrant_count = max(0, int(POP_SIZE * IMMIGRANT_RATE))
        for _ in range(immigrant_count):
            idx = random.randint(0, len(offspring) - 1)
            immigrant = {"codons": mutate(base, protein_seq, mutation_rate=1.0)}
            update_scores(immigrant, protein_seq=protein_seq)
            offspring[idx] = immigrant

        population = select_next_generation_nsga2(population + offspring, POP_SIZE)

#############################################
# Final outputs
#############################################
if MODE == "scalar":
    population.sort(key=lambda x: x["fitness"], reverse=True)
else:
    rank_population(population)
    population.sort(key=lambda x: (x["rank"], -x["crowding"], -x["fitness"]))

best = population[0]
best_seq = best["dna"]

with open(output_fasta, "w") as f:
    f.write(f">{record.id}_optimized\n{best_seq}\n")

# Pareto export (in scalar mode, export top-k by fitness)
if MODE == "pareto":
    pareto_front = [ind for ind in population if ind.get("rank", 9999) == 0]
    pareto_front.sort(key=lambda x: x["fitness"], reverse=True)
else:
    pareto_front = sorted(population, key=lambda x: x["fitness"], reverse=True)

top_k = pareto_front[:PARETO_TOP_K]
pareto_rows = []
for i, ind in enumerate(top_k, start=1):
    pareto_rows.append(
        {
            "rank_in_export": i,
            "fitness": ind["fitness"],
            "cai_r": ind["cai_r"],
            "cai_g": ind["cai_g"],
            "gc": ind["gc"],
            "gc_target": GC_TARGET,
            "gc_score": ind["gc_score"],
            "splice_score": ind["splice_score"],
            "sequence_length": len(ind["dna"]),
            "mode": MODE,
            "pareto_rank": ind.get("rank", 0 if MODE == "scalar" else None),
            "crowding": ind.get("crowding", None),
        }
    )

pd.DataFrame(pareto_rows).to_csv(pareto_front_tsv, sep="\t", index=False)

with open(pareto_sequences_fasta, "w") as f:
    for i, ind in enumerate(top_k, start=1):
        f.write(
            f">{record.id}_{MODE}_candidate_{i}"
            f"_fitness_{ind['fitness']:.6f}"
            f"_caiR_{ind['cai_r']:.4f}"
            f"_caiG_{ind['cai_g']:.4f}"
            f"_gc_{ind['gc']:.4f}"
            "\n"
        )
        f.write(ind["dna"] + "\n")

#############################################
# Write GA Log
#############################################
ga_log_df = pd.DataFrame(log_rows)
ga_log_df["last_generation"] = last_generation
ga_log_df.to_csv(output_log, sep="\t", index=False)

#############################################
# Fitness vs Generation Plot
#############################################
plt.figure()
plt.plot(ga_log_df["generation"], ga_log_df["best_fitness"], label="Best")
plt.plot(ga_log_df["generation"], ga_log_df["mean_fitness"], label="Mean")
if "pareto_front_size" in ga_log_df.columns and MODE == "pareto":
    plt.plot(ga_log_df["generation"], ga_log_df["pareto_front_size"], label="Front size")
plt.xlabel("Generation")
plt.ylabel("Metric")
plt.title(f"GA Progress ({MODE})")
plt.legend()
plt.tight_layout()
plt.savefig(output_plot)
plt.close()