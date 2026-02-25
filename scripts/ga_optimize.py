#!/usr/bin/env python3
"""
Codon optimisation GA — v3
==========================

Snakemake inputs
----------------
REQUIRED:
  snakemake.input.protein         — protein FASTA (single record)
  snakemake.input.ribo_weights    — JSON ribosomal codon weights
  snakemake.input.genome_weights  — JSON genome-wide codon weights
  snakemake.input.genome_cds      — CDS FASTA used for GC target + codon pairs
  snakemake.input.codon_pairs     — JSON codon-pair log-odds table
                                    (build with scripts/build_codon_pair_table.py)

OPTIONAL:
  snakemake.input.homologue_cds   — FASTA of homologous CDS sequences for
                                    population seeding

Config sections (all keys shown with defaults)
----------------------------------------------
ga:
  mode: "pareto"            # scalar | pareto
  pareto_top_k: 10
  population_size: 120
  generations: 10000
  mutation_rate: 0.07
  mutation_rate_final: 0.02
  elite_size: 10
  crossover_rate: 0.75
  tournament_size: 4
  immigrant_rate: 0.05
  early_stop_patience: 1000
  min_improvement: 1.0e-10
  random_seed: 42
  n_islands: 4
  migration_interval: 25
  migration_size: 2
  adaptive_mutation: true
  adapt_window: 50
  n_crossover_points: 2
  max_homologue_seeds: 10

fitness:
  alpha: 0.7                # CAI_ribo exponent
  beta: 0.3                 # CAI_genome exponent
  ramp_length: 50           # codons in 5-prime ramp
  ramp_weight: 0.5          # down-weight factor for ramp codons
  codon_pair_weight: 1.0    # exponent on codon-pair score (0 = disable)

gc:
  use_genome_target: true
  target: 0.52              # fallback GC target
  sigma: 0.05               # Gaussian half-width for GC penalty

splice:
  donor_penalty: 0.8        # per-hit multiplier for splice-donor motifs
  intron_penalty: 0.9       # per-hit multiplier for GT..AG intron-like runs
  intron_min_distance: 20
  intron_max_distance: 200

fungal:
  splice_penalty_weight: 5  # higher = stronger fungal splice-site avoidance
  polyA_penalty_weight: 3   # higher = stronger polyA signal avoidance
  gc_window: 50             # nt window for local GC scan
  low_gc: 0.30              # lower bound of acceptable local GC
  high_gc: 0.70             # upper bound of acceptable local GC
  polya_motifs: ["AATAAA", "ATTAAA"]
  donor_motifs: ["AGGT", "CAGGT", "AAGGT", "GTATGT"]

snapshot:
  count: 10
"""

import json
import math
import random
import warnings
import collections
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
protein_fasta       = snakemake.input.protein
ribo_weights_file   = snakemake.input.ribo_weights
genome_weights_file = snakemake.input.genome_weights
genome_cds_fasta    = snakemake.input.genome_cds
codon_pairs_file    = snakemake.input.codon_pairs
_hom                = getattr(snakemake.input, "homologue_cds", None)
homologue_cds_fasta = _hom if _hom else None

output_fasta           = snakemake.output.optimized
output_log             = snakemake.output.log
output_plot            = snakemake.output.plot
snapshots_fasta        = snakemake.output.snapshots
pareto_front_tsv       = snakemake.output.pareto_front
pareto_sequences_fasta = snakemake.output.pareto_sequences
config                 = snakemake.config

open(snapshots_fasta, "w").close()

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
ga_cfg = config["ga"]
MODE   = ga_cfg.get("mode", "scalar").lower()
if MODE not in {"scalar", "pareto"}:
    raise ValueError("ga.mode must be 'scalar' or 'pareto'")

POP_SIZE            = int(ga_cfg["population_size"])
GENERATIONS         = int(ga_cfg["generations"])
MUTATION_RATE       = float(ga_cfg["mutation_rate"])
ELITE_SIZE          = int(ga_cfg["elite_size"])
RANDOM_SEED         = int(ga_cfg["random_seed"])
PARETO_TOP_K        = max(1, int(ga_cfg.get("pareto_top_k", 10)))
CROSSOVER_RATE      = float(ga_cfg.get("crossover_rate", 0.75))
TOURNAMENT_SIZE     = int(ga_cfg.get("tournament_size", 3))
IMMIGRANT_RATE      = float(ga_cfg.get("immigrant_rate", 0.05))
EARLY_STOP_PATIENCE = int(ga_cfg.get("early_stop_patience", 1000))
MIN_IMPROVEMENT     = float(ga_cfg.get("min_improvement", 1e-10))
MUTATION_RATE_FINAL = float(ga_cfg.get("mutation_rate_final", MUTATION_RATE))

N_ISLANDS          = int(ga_cfg.get("n_islands", 4))
MIGRATION_INTERVAL = int(ga_cfg.get("migration_interval", 25))
MIGRATION_SIZE     = int(ga_cfg.get("migration_size", 2))

USE_ADAPTIVE_MUTATION = bool(ga_cfg.get("adaptive_mutation", True))
ADAPT_WINDOW          = int(ga_cfg.get("adapt_window", 50))
N_CROSSOVER_POINTS    = int(ga_cfg.get("n_crossover_points", 2))
MAX_HOMOLOGUE_SEEDS   = int(ga_cfg.get("max_homologue_seeds", 10))
MUTATION_GAMMA = float(ga_cfg.get("mutation_bias_gamma", 0.0))
random.seed(RANDOM_SEED)

fit_cfg           = config.get("fitness", {})
alpha             = float(fit_cfg.get("alpha", 0.7))
beta              = float(fit_cfg.get("beta", 0.3))
RAMP_LENGTH       = int(fit_cfg.get("ramp_length", 50))
RAMP_WEIGHT       = float(fit_cfg.get("ramp_weight", 0.5))
CODON_PAIR_WEIGHT = float(fit_cfg.get("codon_pair_weight", 1.0))
# Per-term exponent weights — raise each penalty term to this power so its
# influence on the scalar fitness is directly comparable to alpha/beta.
# Defaults to 1.0 (unchanged behaviour). Reduce below 1.0 to soften a term,
# increase above 1.0 to strengthen it relative to CAI.
GC_WEIGHT       = float(fit_cfg.get("gc_weight",       1.0))
SPLICE_WEIGHT   = float(fit_cfg.get("splice_weight",   1.0))
POLYA_WEIGHT    = float(fit_cfg.get("polya_weight",    1.0))
LOCAL_GC_WEIGHT = float(fit_cfg.get("local_gc_weight", 1.0))

gc_cfg               = config.get("gc", {})
USE_GENOME_GC_TARGET = bool(gc_cfg.get("use_genome_target", True))
GC_TARGET_FALLBACK   = float(gc_cfg.get("target", 0.5))
GC_SIGMA             = float(gc_cfg.get("sigma", 0.05))

splice_cfg     = config.get("splice", {})
DONOR_PENALTY  = float(splice_cfg.get("donor_penalty", 0.8))
INTRON_PENALTY = float(splice_cfg.get("intron_penalty", 0.9))
INTRON_MIN     = int(splice_cfg.get("intron_min_distance", 20))
INTRON_MAX     = int(splice_cfg.get("intron_max_distance", 200))

# Fungal-specific parameters — all read from config["fungal"]
fungal_cfg             = config.get("fungal", {})
FUNGAL_SPLICE_WEIGHT   = float(fungal_cfg.get("splice_penalty_weight", 5))
FUNGAL_POLYA_WEIGHT    = float(fungal_cfg.get("polyA_penalty_weight", 3))
FUNGAL_GC_WINDOW       = int(fungal_cfg.get("gc_window", 50))
FUNGAL_LOW_GC          = float(fungal_cfg.get("low_gc", 0.30))
FUNGAL_HIGH_GC         = float(fungal_cfg.get("high_gc", 0.70))
FUNGAL_POLYA_MOTIFS    = list(fungal_cfg.get("polya_motifs", ["AATAAA", "ATTAAA"]))
FUNGAL_DONOR_MOTIFS    = list(fungal_cfg.get("donor_motifs",
                              ["AGGT", "CAGGT", "AAGGT", "GTATGT"]))



snapshot_cfg   = config.get("snapshot", config.get("Snapshot", {}))
SNAPSHOT_COUNT = max(1, int(snapshot_cfg.get("count", 10)))
# Base interval on early_stop_patience so we always get SNAPSHOT_COUNT
# snapshots even when the GA stops well before GENERATIONS.
_snapshot_horizon = min(GENERATIONS, EARLY_STOP_PATIENCE)
SN = max(1, _snapshot_horizon // SNAPSHOT_COUNT)

# ---------------------------------------------------------------------------
# Load weight models
# ---------------------------------------------------------------------------
with open(ribo_weights_file) as f:
    ribo_weights = json.load(f)

with open(genome_weights_file) as f:
    genome_weights = json.load(f)

with open(codon_pairs_file) as f:
    codon_pair_scores = json.load(f)

# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------
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
    "*": ["TAA", "TAG", "TGA"],
}

codon_to_aa = {
    codon: aa
    for aa, codons in aa_to_codons.items()
    for codon in codons
}

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def compute_gc(seq):
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_genome_gc_target(cds_fasta_path, fallback):
    total_gc = total_len = 0
    for rec in SeqIO.parse(cds_fasta_path, "fasta"):
        s = str(rec.seq).upper()
        total_gc  += s.count("G") + s.count("C")
        total_len += len(s)
    return total_gc / total_len if total_len else fallback


GC_TARGET = (
    compute_genome_gc_target(genome_cds_fasta, GC_TARGET_FALLBACK)
    if USE_GENOME_GC_TARGET else GC_TARGET_FALLBACK
)


def backtranslate_max(protein_seq):
    codons, skipped = [], []
    for pos, aa in enumerate(protein_seq):
        if aa not in aa_to_codons:
            skipped.append((pos, aa))
            continue
        codons.append(max(aa_to_codons[aa], key=lambda c: ribo_weights.get(c, 1e-8)))
    if skipped:
        warnings.warn(
            "backtranslate_max: skipped unrecognised residue(s): "
            + ", ".join(f"pos {p} '{a}'" for p, a in skipped)
        )
    return codons


def codon_list_to_seq(codon_list):
    return "".join(codon_list)

# ---------------------------------------------------------------------------
# Synonymy guard
# ---------------------------------------------------------------------------
def verify_synonymous(codon_list, protein_seq, context=""):
    if len(codon_list) != len(protein_seq):
        raise ValueError(
            f"[{context}] Length mismatch: {len(codon_list)} codons vs "
            f"{len(protein_seq)} residues."
        )
    mismatches = []
    for i, (codon, aa) in enumerate(zip(codon_list, protein_seq)):
        if codon_to_aa.get(codon) != aa:
            mismatches.append(
                f"  pos {i}: codon '{codon}' -> '{codon_to_aa.get(codon)}', "
                f"expected '{aa}'"
            )
    if mismatches:
        raise ValueError(
            f"[{context}] Non-synonymous mutation(s) detected "
            f"({len(mismatches)} site(s)):\n"
            + "\n".join(mismatches[:10])
            + ("\n  ..." if len(mismatches) > 10 else "")
        )

# ---------------------------------------------------------------------------
# Fitness components
# ---------------------------------------------------------------------------
def compute_cai_windowed(codon_list, weights):
    """CAI with a 5-prime translational ramp."""
    log_sum = eff_len = 0.0
    for i, codon in enumerate(codon_list):
        w          = max(0.01, (weights.get(codon, 1e-8) - 0.1))
        pos_weight = RAMP_WEIGHT if i < RAMP_LENGTH else 1.0
        log_sum   += pos_weight * math.log(w)
        eff_len   += pos_weight
    return math.exp(log_sum / eff_len) if eff_len else 0.0


def compute_codon_pair_score(codon_list):
    """Sigmoid-normalised mean codon-pair log-odds score."""
    if len(codon_list) < 2 or CODON_PAIR_WEIGHT == 0:
        return 1.0
    total = sum(
        codon_pair_scores.get(codon_list[i] + codon_list[i + 1], 0.0)
        for i in range(len(codon_list) - 1)
    )
    mean_log_odds = total / (len(codon_list) - 1)
    return 1.0 / (1.0 + math.exp(-mean_log_odds))


def gc_penalty(seq):
    """Global GC Gaussian penalty centred on GC_TARGET."""
    gc = compute_gc(seq)
    return math.exp(-((gc - GC_TARGET) ** 2) / (2 * GC_SIGMA ** 2))


def fungal_splice_penalty(seq):
    """
    Fungal splice-site avoidance.

    Uses the motifs and weights from config['fungal'] rather than the
    generic splice block.  Each motif hit multiplies the score by
    DONOR_PENALTY raised to (1 / FUNGAL_SPLICE_WEIGHT), so higher
    FUNGAL_SPLICE_WEIGHT = steeper penalty per hit.

    Additionally checks for GT..AG pseudo-intron windows as before.
    The per-hit penalty is scaled by FUNGAL_SPLICE_WEIGHT so that the
    config value directly controls how aggressively these sites are avoided.
    """
    penalty = 1.0
    length  = len(seq)

    # Motif-based donor penalties — uses fungal donor_motifs from config
    per_hit = DONOR_PENALTY ** (1.0 / max(FUNGAL_SPLICE_WEIGHT, 1e-6))
    for motif in FUNGAL_DONOR_MOTIFS:
        count = seq.count(motif)
        if count:
            penalty *= per_hit ** count

    # GT..AG intron-like window scan
    per_intron = INTRON_PENALTY ** (1.0 / max(FUNGAL_SPLICE_WEIGHT, 1e-6))
    for i in range(length - 2):
        if seq[i: i + 2] == "GT":
            window = seq[i + INTRON_MIN: i + INTRON_MAX]
            if "AG" in window:
                penalty *= per_intron

    return penalty


def fungal_polya_penalty(seq):
    """
    Penalise cryptic polyadenylation signals (AATAAA, ATTAAA and any
    additional motifs from config['fungal']['polya_motifs']).

    Each hit multiplies the score by exp(-1 / FUNGAL_POLYA_WEIGHT), so a
    higher weight = stronger avoidance.  The exponential form keeps the
    penalty in (0, 1] regardless of hit count.
    """
    penalty   = 1.0
    per_hit   = math.exp(-1.0 / max(FUNGAL_POLYA_WEIGHT, 1e-6))
    for motif in FUNGAL_POLYA_MOTIFS:
        count = seq.count(motif)
        if count:
            penalty *= per_hit ** count
    return penalty

def fungal_local_gc_penalty(seq):
    """
    Sliding window GC penalty (AVERAGED, not multiplicative).

    Returns mean window penalty in (0,1].
    Prevents exponential collapse for long CDS.
    """
    w = FUNGAL_GC_WINDOW
    length = len(seq)

    if length < w:
        return 1.0

    penalties = []

    # Pre-count GC in first window
    gc_count = seq[:w].count("G") + seq[:w].count("C")

    for i in range(length - w + 1):

        if i > 0:
            out_base = seq[i - 1]
            in_base = seq[i + w - 1]

            if out_base in "GC":
                gc_count -= 1
            if in_base in "GC":
                gc_count += 1

        local_gc = gc_count / w

        if local_gc < FUNGAL_LOW_GC:
            delta = FUNGAL_LOW_GC - local_gc
            penalties.append(
                math.exp(-(delta ** 2) / (2 * GC_SIGMA ** 2))
            )

        elif local_gc > FUNGAL_HIGH_GC:
            delta = local_gc - FUNGAL_HIGH_GC
            penalties.append(
                math.exp(-(delta ** 2) / (2 * GC_SIGMA ** 2))
            )
        else:
            penalties.append(1.0)

    return sum(penalties) / len(penalties)


# ---------------------------------------------------------------------------
# Score cache + evaluate
# ---------------------------------------------------------------------------

score_cache: dict = {}

def evaluate(codon_list):
    dna = codon_list_to_seq(codon_list)

    cached = score_cache.get(dna)
    if cached is not None:
        return cached

    cai_r       = compute_cai_windowed(codon_list, ribo_weights)
    cai_g       = compute_cai_windowed(codon_list, genome_weights)
    gc_score    = gc_penalty(dna)
    splice_sc   = fungal_splice_penalty(dna)
    polya_sc    = fungal_polya_penalty(dna)
    local_gc_sc = fungal_local_gc_penalty(dna)
    cp_score    = compute_codon_pair_score(codon_list)

    EPS = 1e-12  # prevents log(0)

    # ------------------------------
    # LOG-DOMAIN SCALAR FITNESS
    # ------------------------------
    log_fitness = (
        alpha * math.log(cai_r + EPS)
        + beta * math.log(cai_g + EPS)
        + GC_WEIGHT * math.log(gc_score + EPS)
        + SPLICE_WEIGHT * math.log(splice_sc + EPS)
        + POLYA_WEIGHT * math.log(polya_sc + EPS)
        + LOCAL_GC_WEIGHT * math.log(local_gc_sc + EPS)
        + CODON_PAIR_WEIGHT * math.log(cp_score + EPS)
    )

    result = {
        "dna": dna,
        "cai_r": cai_r,
        "cai_g": cai_g,
        "gc": compute_gc(dna),
        "gc_score": gc_score,
        "splice_score": splice_sc,
        "polya_score": polya_sc,
        "local_gc_score": local_gc_sc,
        "cp_score": cp_score,
        "fitness": log_fitness,  # ← NOW LOG FITNESS
    }

    score_cache[dna] = result
    return result
# ---------------------------------------------------------------------------
# Mutation
# ---------------------------------------------------------------------------
def mutate(codon_list, protein_seq, mutation_rate, position_rates=None):
    new = codon_list.copy()
    for i, aa in enumerate(protein_seq):
        if aa not in aa_to_codons:
            continue
        rate = position_rates[i] if position_rates is not None else mutation_rate
        if random.random() < rate:
            codons = aa_to_codons[aa]

            if MUTATION_GAMMA > 0:
                weights = [
                    ribo_weights.get(c, 1e-8) ** MUTATION_GAMMA
                    for c in codons
                ]
                new[i] = random.choices(codons, weights=weights, k=1)[0]
            else:
                new[i] = random.choice(codons)
    return new


def update_scores(ind, protein_seq=None, debug=True):
    if debug and protein_seq is not None:
        verify_synonymous(ind["codons"], protein_seq, context="update_scores")
    ind.update(evaluate(ind["codons"]))


class AdaptiveMutationTracker:
    def __init__(self, n_positions, base_rate, window=50):
        self.n        = n_positions
        self.base     = base_rate
        self._history = collections.deque(maxlen=window)

    def update(self, population):
        snapshot = []
        for pos in range(self.n):
            counts = collections.Counter(
                ind["codons"][pos] for ind in population
                if pos < len(ind["codons"])
            )
            total = sum(counts.values())
            if total == 0:
                snapshot.append(0.0)
            else:
                top_freq = counts.most_common(1)[0][1] / total
                snapshot.append(1.0 - top_freq)
        self._history.append(snapshot)

    def rates(self):
        if not self._history:
            return None
        n_gen = len(self._history)
        return [
            min(self.base * (0.5 + sum(
                self._history[g][pos] for g in range(n_gen)
            ) / n_gen), 1.0)
            for pos in range(self.n)
        ]

# ---------------------------------------------------------------------------
# Crossover — multi-point
# ---------------------------------------------------------------------------
def crossover(parent_a, parent_b, n_points=None):
    n_points = N_CROSSOVER_POINTS if n_points is None else n_points
    L = len(parent_a)
    if L <= 1:
        return parent_a.copy(), parent_b.copy()
    n_points = min(n_points, L - 1)
    cuts     = sorted(random.sample(range(1, L), n_points))
    child_a, child_b = [], []
    prev, use_ab = 0, True
    for cut in cuts + [L]:
        sa, sb = parent_a[prev:cut], parent_b[prev:cut]
        if use_ab:
            child_a.extend(sa); child_b.extend(sb)
        else:
            child_a.extend(sb); child_b.extend(sa)
        use_ab = not use_ab
        prev   = cut
    return child_a, child_b

# ---------------------------------------------------------------------------
# Scheduling
# ---------------------------------------------------------------------------
def scheduled_mutation_rate(gen):
    if GENERATIONS <= 1:
        return MUTATION_RATE
    progress = gen / (GENERATIONS - 1)
    return MUTATION_RATE + (MUTATION_RATE_FINAL - MUTATION_RATE) * progress

# ---------------------------------------------------------------------------
# NSGA-II
# ---------------------------------------------------------------------------
def dominates(a, b):
    # cp_score added as a 5th objective
    obj_a = (a["cai_r"], a["cai_g"], a["gc_score"], a["splice_score"],
             a["polya_score"], a["local_gc_score"], a["cp_score"])
    obj_b = (b["cai_r"], b["cai_g"], b["gc_score"], b["splice_score"],
             b["polya_score"], b["local_gc_score"], b["cp_score"])
    return (all(x >= y for x, y in zip(obj_a, obj_b))
            and any(x > y for x, y in zip(obj_a, obj_b)))


def fast_non_dominated_sort(pop):
    for p in pop:
        p["dominated_set"] = []
        p["dom_count"]     = 0
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
        nxt = []
        for p in fronts[i]:
            for q in p["dominated_set"]:
                q["dom_count"] -= 1
                if q["dom_count"] == 0:
                    q["rank"] = i + 1
                    nxt.append(q)
        i += 1
        if nxt:
            fronts.append(nxt)
    return fronts


def assign_crowding_distance(front):
    if not front:
        return
    for ind in front:
        ind["crowding"] = 0.0
    objectives = ["cai_r", "cai_g", "gc_score", "splice_score",
                  "polya_score", "local_gc_score", "cp_score"]
    n = len(front)
    if n <= 2:
        for ind in front:
            ind["crowding"] = float("inf")
        return
    for obj in objectives:
        front.sort(key=lambda x: x[obj])
        front[0]["crowding"] = float("inf")
        front[-1]["crowding"] = float("inf")
        obj_min, obj_max = front[0][obj], front[-1][obj]
        if obj_max == obj_min:
            continue
        for i in range(1, n - 1):
            if math.isinf(front[i]["crowding"]):
                continue
            front[i]["crowding"] += (
                (front[i + 1][obj] - front[i - 1][obj]) / (obj_max - obj_min)
            )


def rank_population(pop):
    fronts = fast_non_dominated_sort(pop)
    for front in fronts:
        assign_crowding_distance(front)
    return fronts


def crowded_tournament(pop, k=2):
    contestants = random.sample(pop, min(len(pop), max(2, k)))
    contestants.sort(
        key=lambda x: (x.get("rank", 10**9), -x.get("crowding", 0.0))
    )
    return contestants[0]


def scalar_tournament(pop, k=3):
    k = max(2, min(k, len(pop)))
    return max(random.sample(pop, k), key=lambda ind: ind["fitness"])


def select_next_generation_nsga2(combined, target_size):
    fronts = rank_population(combined)
    new_pop = []
    for front in fronts:
        if len(new_pop) + len(front) <= target_size:
            new_pop.extend(front)
        else:
            front.sort(key=lambda x: x.get("crowding", 0.0), reverse=True)
            new_pop.extend(front[: target_size - len(new_pop)])
            break
    return new_pop

# ---------------------------------------------------------------------------
# Homologue seeding
# ---------------------------------------------------------------------------
def load_homologue_seeds(fasta_path, protein_seq, max_seeds):
    seeds = []
    if not fasta_path:
        return seeds
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if len(seeds) >= max_seeds:
            break
        dna = str(rec.seq).upper().replace("U", "T")
        if len(dna) % 3 != 0:
            continue
        codon_list = [dna[i: i + 3] for i in range(0, len(dna), 3)]
        if codon_list and codon_to_aa.get(codon_list[-1]) == "*":
            codon_list = codon_list[:-1]
        if len(codon_list) != len(protein_seq):
            continue
        if all(codon_to_aa.get(c) == aa
               for c, aa in zip(codon_list, protein_seq)):
            seeds.append(codon_list)
    return seeds

# ---------------------------------------------------------------------------
# Island model
# ---------------------------------------------------------------------------
def migrate(islands, migration_size, mode):
    n         = len(islands)
    emigrants = []
    for island in islands:
        if mode == "scalar":
            island.sort(key=lambda x: x["fitness"], reverse=True)
        else:
            rank_population(island)
            island.sort(
                key=lambda x: (x.get("rank", 9999), -x.get("crowding", 0.0))
            )
        emigrants.append([dict(ind) for ind in island[:migration_size]])

    for i, island in enumerate(islands):
        arrivals = emigrants[(i - 1) % n]
        if mode == "scalar":
            island.sort(key=lambda x: x["fitness"])
        else:
            island.sort(
                key=lambda x: (-x.get("rank", 0), x.get("crowding", 0.0))
            )
        for j, newcomer in enumerate(arrivals):
            island[j] = newcomer


def evolve_one_generation(population, protein_seq, mutation_rate,
                           position_rates, mode):
    island_size = len(population)

    if mode == "scalar":
        new_pop = population[:ELITE_SIZE]
        while len(new_pop) < island_size:
            p1 = scalar_tournament(population, k=TOURNAMENT_SIZE)
            p2 = scalar_tournament(population, k=TOURNAMENT_SIZE)
            c1_cod, c2_cod = (
                crossover(p1["codons"], p2["codons"])
                if random.random() < CROSSOVER_RATE
                else (p1["codons"].copy(), p2["codons"].copy())
            )
            for cod in (c1_cod, c2_cod):
                if len(new_pop) >= island_size:
                    break
                cod   = mutate(cod, protein_seq, mutation_rate, position_rates)
                child = {"codons": cod}
                update_scores(child, protein_seq=protein_seq)
                new_pop.append(child)

        immigrant_count = max(0, int(island_size * IMMIGRANT_RATE))
        for _ in range(immigrant_count):
            if len(new_pop) <= ELITE_SIZE:
                break
            idx = random.randint(ELITE_SIZE, len(new_pop) - 1)
            imm = {"codons": mutate(base, protein_seq, mutation_rate=1.0)}
            update_scores(imm, protein_seq=protein_seq)
            new_pop[idx] = imm
        return new_pop

    else:
        rank_population(population)
        offspring = []
        while len(offspring) < island_size:
            p1 = crowded_tournament(population, k=TOURNAMENT_SIZE)
            p2 = crowded_tournament(population, k=TOURNAMENT_SIZE)
            c1_cod, c2_cod = (
                crossover(p1["codons"], p2["codons"])
                if random.random() < CROSSOVER_RATE
                else (p1["codons"].copy(), p2["codons"].copy())
            )
            for cod in (c1_cod, c2_cod):
                if len(offspring) >= island_size:
                    break
                cod   = mutate(cod, protein_seq, mutation_rate, position_rates)
                child = {"codons": cod}
                update_scores(child, protein_seq=protein_seq)
                offspring.append(child)

        immigrant_count = max(0, int(island_size * IMMIGRANT_RATE))
        for _ in range(immigrant_count):
            idx = random.randint(0, len(offspring) - 1)
            imm = {"codons": mutate(base, protein_seq, mutation_rate=1.0)}
            update_scores(imm, protein_seq=protein_seq)
            offspring[idx] = imm

        return select_next_generation_nsga2(population + offspring, island_size)

# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------
record      = next(SeqIO.parse(protein_fasta, "fasta"))
protein_seq = str(record.seq)
base        = backtranslate_max(protein_seq)
verify_synonymous(base, protein_seq, context="backtranslate_max")

n_codons         = len(base)
island_size      = max(1, POP_SIZE // N_ISLANDS)
homologue_seeds  = load_homologue_seeds(
    homologue_cds_fasta, protein_seq, MAX_HOMOLOGUE_SEEDS
)

def build_initial_island(seed_pool):
    island = []
    for codons in seed_pool[:island_size]:
        ind = {"codons": list(codons)}
        update_scores(ind, protein_seq=protein_seq)
        island.append(ind)
    while len(island) < island_size:
        ind = {"codons": mutate(base, protein_seq, mutation_rate=MUTATION_RATE)}
        update_scores(ind, protein_seq=protein_seq)
        island.append(ind)
    return island

seeds_per_island = [[] for _ in range(N_ISLANDS)]
for idx, seed in enumerate(homologue_seeds):
    seeds_per_island[idx % N_ISLANDS].append(seed)

islands = [build_initial_island(seeds_per_island[i]) for i in range(N_ISLANDS)]

if MODE == "pareto":
    for island in islands:
        rank_population(island)

trackers = [
    AdaptiveMutationTracker(n_codons, MUTATION_RATE, window=ADAPT_WINDOW)
    for _ in range(N_ISLANDS)
]

# ---------------------------------------------------------------------------
# GA loop
# ---------------------------------------------------------------------------
log_rows             = []
best_overall_fitness = -1.0
stagnation_counter   = 0
last_generation      = GENERATIONS - 1

for gen in range(GENERATIONS):
    mutation_rate = scheduled_mutation_rate(gen)

    for i, island in enumerate(islands):
        if USE_ADAPTIVE_MUTATION:
            trackers[i].update(island)
            pos_rates = trackers[i].rates()
        else:
            pos_rates = None
        islands[i] = evolve_one_generation(
            island, protein_seq, mutation_rate, pos_rates, MODE
        )

    if N_ISLANDS > 1 and gen > 0 and gen % MIGRATION_INTERVAL == 0:
        migrate(islands, MIGRATION_SIZE, MODE)

    all_inds = [ind for island in islands for ind in island]
    if MODE == "scalar":
        all_inds.sort(key=lambda x: x["fitness"], reverse=True)
    else:
        rank_population(all_inds)
        all_inds.sort(
            key=lambda x: (x.get("rank", 9999), -x.get("crowding", 0.0))
        )

    best_fit   = all_inds[0]["fitness"]
    mean_fit   = sum(ind["fitness"] for ind in all_inds) / len(all_inds)
    front_size = (
        sum(1 for ind in all_inds if ind.get("rank", 9999) == 0)
        if MODE == "pareto" else 1
    )

    log_rows.append({
        "generation":        gen,
        "best_fitness":      best_fit,
        "mean_fitness":      mean_fit,
        "mutation_rate":     mutation_rate,
        "cache_size":        len(score_cache),
        "mode":              MODE,
        "pareto_front_size": front_size,
        "gc_target":         GC_TARGET,
        "n_islands":         N_ISLANDS,
        # Per-component scores of the current best individual
        "best_cai_r":        all_inds[0].get("cai_r",         0.0),
        "best_cai_g":        all_inds[0].get("cai_g",         0.0),
        "best_gc_score":     all_inds[0].get("gc_score",      0.0),
        "best_splice_score": all_inds[0].get("splice_score",  0.0),
        "best_polya_score":  all_inds[0].get("polya_score",   0.0),
        "best_local_gc":     all_inds[0].get("local_gc_score",0.0),
        "best_cp_score":     all_inds[0].get("cp_score",      0.0),
    })

    if best_fit > best_overall_fitness + MIN_IMPROVEMENT:
        best_overall_fitness = best_fit
        stagnation_counter   = 0
    else:
        stagnation_counter += 1

    if gen % SN == 0 or gen == GENERATIONS - 1:
        best_seq = all_inds[0]["dna"]
        with open(snapshots_fasta, "a") as f:
            f.write(f">generation_{gen}_fitness_{best_fit:.6f}\n{best_seq}\n")

    if stagnation_counter >= EARLY_STOP_PATIENCE:
        last_generation = gen
        break

# ---------------------------------------------------------------------------
# Final outputs
# ---------------------------------------------------------------------------
all_inds = [ind for island in islands for ind in island]

if MODE == "scalar":
    all_inds.sort(key=lambda x: x["fitness"], reverse=True)
else:
    rank_population(all_inds)
    all_inds.sort(
        key=lambda x: (x.get("rank", 9999), -x.get("crowding", 0.0), -x["fitness"])
    )

best     = all_inds[0]
best_seq = best["dna"]

with open(output_fasta, "w") as f:
    f.write(f">{record.id}_optimized\n{best_seq}\n")

pareto_front = (
    sorted(
        [ind for ind in all_inds if ind.get("rank", 9999) == 0],
        key=lambda x: x["fitness"], reverse=True
    )
    if MODE == "pareto"
    else sorted(all_inds, key=lambda x: x["fitness"], reverse=True)
)

top_k       = pareto_front[:PARETO_TOP_K]
pareto_rows = []
for i, ind in enumerate(top_k, start=1):
    pareto_rows.append({
        "rank_in_export":   i,
        "fitness":          ind["fitness"],
        "cai_r":            ind["cai_r"],
        "cai_g":            ind["cai_g"],
        "gc":               ind["gc"],
        "gc_target":        GC_TARGET,
        "gc_score":         ind["gc_score"],
        "splice_score":     ind["splice_score"],
        "polya_score":      ind["polya_score"],
        "local_gc_score":   ind["local_gc_score"],
        "cp_score":         ind["cp_score"],
        "sequence_length":  len(ind["dna"]),
        "mode":             MODE,
        "pareto_rank":      ind.get("rank", 0 if MODE == "scalar" else None),
        "crowding":         ind.get("crowding", None),
    })

pd.DataFrame(pareto_rows).to_csv(pareto_front_tsv, sep="\t", index=False)

with open(pareto_sequences_fasta, "w") as f:
    for i, ind in enumerate(top_k, start=1):
        f.write(
            f">{record.id}_{MODE}_candidate_{i}"
            f"_fitness_{ind['fitness']:.6f}"
            f"_caiR_{ind['cai_r']:.4f}"
            f"_caiG_{ind['cai_g']:.4f}"
            f"_gc_{ind['gc']:.4f}"
            f"_polyA_{ind['polya_score']:.4f}"
            f"_localGC_{ind['local_gc_score']:.4f}"
            f"_cp_{ind['cp_score']:.4f}"
            "\n"
        )
        f.write(ind["dna"] + "\n")

# ---------------------------------------------------------------------------
# Log
# ---------------------------------------------------------------------------
ga_log_df = pd.DataFrame(log_rows)
ga_log_df["last_generation"] = last_generation
ga_log_df.to_csv(output_log, sep="\t", index=False)

# ---------------------------------------------------------------------------
# Plot — main fitness + per-component breakdown
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax = axes[0]
ax.plot(ga_log_df["generation"], ga_log_df["best_fitness"], label="Best fitness")
ax.plot(ga_log_df["generation"], ga_log_df["mean_fitness"], label="Mean fitness",
        alpha=0.6)
if MODE == "pareto" and "pareto_front_size" in ga_log_df.columns:
    ax2 = ax.twinx()
    ax2.plot(ga_log_df["generation"], ga_log_df["pareto_front_size"],
             color="grey", alpha=0.4, linestyle="--", label="Front size")
    ax2.set_ylabel("Pareto front size")
ax.set_ylabel("Fitness")
ax.set_title(f"GA Progress ({MODE}, {N_ISLANDS} islands)")
ax.legend(loc="lower right")

ax = axes[1]
for col, label in [
    ("best_cai_r",        "CAI ribo"),
    ("best_cai_g",        "CAI genome"),
    ("best_gc_score",     "Global GC"),
    ("best_splice_score", "Splice"),
    ("best_polya_score",  "PolyA"),
    ("best_local_gc",     "Local GC"),
    ("best_cp_score",     "Codon-pair"),
]:
    if col in ga_log_df.columns:
        ax.plot(ga_log_df["generation"], ga_log_df[col], label=label)
ax.set_xlabel("Generation")
ax.set_ylabel("Component score")
ax.set_title("Best individual — per-component scores")
ax.legend(loc="lower right", fontsize=8)

plt.tight_layout()
plt.savefig(output_plot)
plt.close()
