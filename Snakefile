#############################################
# Codon Optimization Pipeline
# Organism: Aureobasidium pullulans
#############################################

configfile: "config.yaml"

import os
from datetime import datetime
import sys

#############################################
# RUN ID LOGIC
#############################################

RUN_ID = config.get("run", {}).get("id")

if not RUN_ID:
    RUN_ID = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

RESULTS_DIR = f"results/{RUN_ID}"
sys.stderr.write(f"RUN_ID: {RUN_ID}\n")
sys.stderr.write(f"Results directory: {RESULTS_DIR}\n")



#############################################
# FINAL TARGET
#############################################

rule all:
    input:
        f"{RESULTS_DIR}/{RUN_ID}optimization_report.pdf"

#############################################
# -------------------------------
# PHASE 1: GENOME CODON ANALYSIS
# -------------------------------
#############################################

rule save_config_snapshot:
    output:
        f"{RESULTS_DIR}/config_used.yaml"
    run:
        import shutil
        shutil.copy("config.yaml", output[0])

rule extract_cds:
    input:
        config_snapshot=f"{RESULTS_DIR}/config_used.yaml",
        genome=config["genome"],
        gff=config["gff"]
    conda:
        "Reg"
    output:
        f"{RESULTS_DIR}/all_cds.fasta"
    script:
        "scripts/extract_cds.py"

rule codon_bias_analysis:
    input:
        f"{RESULTS_DIR}/all_cds.fasta"
    conda:
        "bio-r"
    output:
        metrics=f"{RESULTS_DIR}/gene_metrics.tsv",
        wright=f"{RESULTS_DIR}/wright_plot.png"
    script:
        "scripts/codon_bias_analysis.R"

rule identify_ribosomal:
    input:
        gff=config["gff"]
    output:
        f"{RESULTS_DIR}/ribosomal_ids.txt"
    script:
        "scripts/find_ribosomal.py"

rule build_ribosomal_model:
    input:
        cds=f"{RESULTS_DIR}/all_cds.fasta",
        ids=f"{RESULTS_DIR}/ribosomal_ids.txt"
    conda:
        "Reg"
    output:
        f"{RESULTS_DIR}/codon_weights_ribo.json"
    script:
        "scripts/build_cai_model.py"

rule build_genome_model:
    input:
        cds=f"{RESULTS_DIR}/all_cds.fasta"
    conda:
        "Reg"
    output:
        f"{RESULTS_DIR}/codon_weights_genome.json"
    script:
        "scripts/build_cai_model.py"
#############################################
# -------------------------------
# PHASE 2: CONSTRUCT OPTIMIZATION
# -------------------------------
#############################################
rule ga_optimize:
    input:
        protein=config["protein"],
        ribo_weights=f"{RESULTS_DIR}/codon_weights_ribo.json",
        genome_weights=f"{RESULTS_DIR}/codon_weights_genome.json",
        genome_cds=f"{RESULTS_DIR}/all_cds.fasta"
    conda:
        "Reg"
    output:
        optimized=f"{RESULTS_DIR}/optimized_sequence.fasta",
        log=f"{RESULTS_DIR}/ga_log.tsv",
        plot=f"{RESULTS_DIR}/fitness_plot.png",
        snapshots=f"{RESULTS_DIR}/intermediary_snapshots.fasta",
        pareto_front=f"{RESULTS_DIR}/pareto_front.tsv",
        pareto_sequences=f"{RESULTS_DIR}/pareto_sequences.fasta"
    script:
        "scripts/ga_optimize.py"

rule fungal_filter:
    input:
        f"{RESULTS_DIR}/optimized_sequence.fasta"
    conda:
        "Reg"
    output:
        f"{RESULTS_DIR}/final_sequence.fasta",
        f"{RESULTS_DIR}/fungal_warnings.txt"
    script:
        "scripts/fungal_filters.py"

rule cai_snapshot_plots:
    input:
        snapshots=f"{RESULTS_DIR}/intermediary_snapshots.fasta",
        genome_weights=f"{RESULTS_DIR}/codon_weights_genome.json"
    output:
        directory(f"{RESULTS_DIR}/cai_plots")
    conda:
        "Reg"
    script:
        "scripts/cai_snapshot_plots.py"


rule codon_usage_grid:
    input:
        genome_weights=f"{RESULTS_DIR}/codon_weights_genome.json",
        ribo_weights=f"{RESULTS_DIR}/codon_weights_ribo.json",
        final_seq=f"{RESULTS_DIR}/final_sequence.fasta"
    output:
        directory(f"{RESULTS_DIR}/codon_usage_grid")
    conda:
        "Reg"
    script:
        "scripts/codon_usage_grid.py"


rule diagnostics_plots:
    input:
        all_cds=f"{RESULTS_DIR}/all_cds.fasta",
        gene_metrics=f"{RESULTS_DIR}/gene_metrics.tsv",
        ribosomal_ids=f"{RESULTS_DIR}/ribosomal_ids.txt",
        ribo_weights=f"{RESULTS_DIR}/codon_weights_ribo.json",
        genome_weights=f"{RESULTS_DIR}/codon_weights_genome.json",
        final_seq=f"{RESULTS_DIR}/final_sequence.fasta",
        pareto_front=f"{RESULTS_DIR}/pareto_front.tsv",
        pareto_sequences=f"{RESULTS_DIR}/pareto_sequences.fasta",
        ga_log=f"{RESULTS_DIR}/ga_log.tsv"
    output:
        directory(f"{RESULTS_DIR}/diagnostics")
    conda:
        "Reg"
    script:
        "scripts/diagnostic_plots.py"

rule final_report:
    input:
        ga_log=f"{RESULTS_DIR}/ga_log.tsv",
        final_seq=f"{RESULTS_DIR}/final_sequence.fasta",
        wright=f"{RESULTS_DIR}/wright_plot.png",
        fitness_plot=f"{RESULTS_DIR}/fitness_plot.png",
        pareto_front=f"{RESULTS_DIR}/pareto_front.tsv",
        diagnostics=directory(f"{RESULTS_DIR}/diagnostics"),
        cai_dir = directory(f"{RESULTS_DIR}/cai_plots"),
        codon_usage_grid = directory(f"{RESULTS_DIR}/codon_usage_grid")
    conda:
        "Reg"
    output:
        f"{RESULTS_DIR}/{RUN_ID}optimization_report.pdf"
    script:
        "scripts/report.py"