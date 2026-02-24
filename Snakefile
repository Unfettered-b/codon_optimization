#############################################
# Codon Optimization Pipeline
# Organism: Aureobasidium pullulans
#############################################

configfile: "config.yaml"

#############################################
# FINAL TARGET
#############################################

rule all:
    input:
        "results/optimization_report.pdf"

#############################################
# -------------------------------
# PHASE 1: GENOME CODON ANALYSIS
# -------------------------------
#############################################

rule extract_cds:
    input:
        genome=config["genome"],
        gff=config["gff"]
    conda:
        "Reg"
    output:
        "results/all_cds.fasta"
    script:
        "scripts/extract_cds.py"

rule codon_bias_analysis:
    input:
        "results/all_cds.fasta"
    conda:
        "bio-r"
    output:
        metrics="results/gene_metrics.tsv",
        wright="results/wright_plot.pdf"
    script:
        "scripts/codon_bias_analysis.R"

rule identify_ribosomal:
    input:
        gff=config["gff"]
    output:
        "results/ribosomal_ids.txt"
    script:
        "scripts/find_ribosomal.py"


rule build_ribosomal_model:
    input:
        cds="results/all_cds.fasta",
        ids="results/ribosomal_ids.txt"
    conda:
        "Reg"
    output:
        "results/codon_weights_ribo.json"
    script:
        "scripts/build_cai_model.py"

rule build_genome_model:
    input:
        cds="results/all_cds.fasta"
    conda:
        "Reg"
    output:
        "results/codon_weights_genome.json"
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
        ribo_weights="results/codon_weights_ribo.json",
        genome_weights="results/codon_weights_genome.json"
    conda:
        "Reg"
    output:
        "results/optimized_sequence.fasta",
        "results/ga_log.tsv"
    script:
        "scripts/ga_optimize.py"


rule fungal_filter:
    input:
        "results/optimized_sequence.fasta"
    conda:
        "Reg"
    output:
        "results/final_sequence.fasta",
        "results/fungal_warnings.txt"
    script:
        "scripts/fungal_filters.py"


rule final_report:
    input:
        ga_log="results/ga_log.tsv",
        final_seq="results/final_sequence.fasta",
        wright="results/wright_plot.pdf"
    conda:
        "Reg"
    output:
        "results/optimization_report.pdf"
    script:
        "scripts/report.py"