#!/usr/bin/env python3

import os
import re

import pandas as pd
from Bio import SeqIO
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
)

#############################################
# Inputs
#############################################

ga_log_file = snakemake.input.ga_log
final_seq_file = snakemake.input.final_seq
wright_plot_file = snakemake.input.wright
fitness_plot_file = snakemake.input.fitness_plot
pareto_front_file = snakemake.input.pareto_front
diagnostics_dir = snakemake.input.diagnostics
cai_dir = snakemake.input.cai_dir
grid_dir = snakemake.input.codon_usage_grid

output_pdf = snakemake.output[0]

#############################################
# Load Data
#############################################

ga_log_df = pd.read_csv(ga_log_file, sep="\t")
pareto_df = pd.read_csv(pareto_front_file, sep="\t")
record = next(SeqIO.parse(final_seq_file, "fasta"))
final_seq = str(record.seq)


def format_fasta(seq):
    lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    return "\n".join(lines)


def extract_generation(filename):
    match = re.search(r"generation_(\d+)", filename)
    if match:
        return int(match.group(1))
    return -1


def compute_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_gc3(seq):
    third = seq[2::3]
    return (third.count("G") + third.count("C")) / len(third)


final_gc = compute_gc(final_seq)
final_gc3 = compute_gc3(final_seq)
best_fitness = ga_log_df["best_fitness"].max()
final_generation = ga_log_df["generation"].max()

#############################################
# PDF Setup
#############################################

doc = SimpleDocTemplate(output_pdf)
elements = []

styles = getSampleStyleSheet()
title_style = styles["Heading1"]
normal_style = styles["Normal"]
heading_style = styles["Heading2"]
sub_heading_style = styles["Heading3"]
fasta_style = ParagraphStyle(name="FastaStyle", fontName="Courier", fontSize=8, leading=10)

#############################################
# Title
#############################################

elements.append(Paragraph("Codon Optimization Report", title_style))
elements.append(Spacer(1, 0.3 * inch))

#############################################
# GA Summary
#############################################

elements.append(Paragraph("1. Genetic Algorithm Optimization Summary", heading_style))
elements.append(Spacer(1, 0.2 * inch))
elements.append(Paragraph(f"Total generations: {final_generation}", normal_style))
elements.append(Paragraph(f"Best composite fitness achieved: {round(best_fitness, 6)}", normal_style))

if "mode" in ga_log_df.columns:
    mode = str(ga_log_df["mode"].iloc[-1])
    elements.append(Paragraph(f"Optimization mode: {mode}", normal_style))

if len(pareto_df) > 0:
    top = pareto_df.iloc[0]
    elements.append(Paragraph(f"Top exported candidates: {len(pareto_df)}", normal_style))
    elements.append(
        Paragraph(
            f"Top candidate metrics — CAI_r: {top['cai_r']:.4f}, CAI_g: {top['cai_g']:.4f}, GC: {top['gc']:.4f}",
            normal_style,
        )
    )

if "gc_target" in ga_log_df.columns:
    elements.append(Paragraph(f"GC target used in optimization: {float(ga_log_df['gc_target'].iloc[-1]):.4f}", normal_style))

elements.append(Spacer(1, 0.3 * inch))

elements.append(Paragraph("1.1 Fitness Progression", sub_heading_style))
elements.append(Spacer(1, 0.1 * inch))
elements.append(Image(fitness_plot_file, width=5 * inch, height=3 * inch))

#############################################
# CAI Snapshot Plots
#############################################

elements.append(PageBreak())
elements.append(Paragraph("1.2 CAI Positional Plots Across Generations", sub_heading_style))
elements.append(Spacer(1, 0.2 * inch))

cai_images = [f for f in os.listdir(cai_dir) if f.endswith(".png")]
cai_images.sort(key=extract_generation)

for img_file in cai_images:
    gen_number = extract_generation(img_file)
    img_path = os.path.join(cai_dir, img_file)
    elements.append(Paragraph(f"Generation {gen_number}", styles["Heading3"]))
    elements.append(Image(img_path, width=6 * inch, height=1.5 * inch))
    elements.append(Spacer(1, 0.15 * inch))

#############################################
# Final sequence + codon grid
#############################################

elements.append(PageBreak())
elements.append(Paragraph("2. Final Optimized Sequence Metrics", heading_style))
elements.append(Spacer(1, 0.2 * inch))
elements.append(Paragraph(f"Sequence length: {len(final_seq)} nt", normal_style))
elements.append(Paragraph(f"GC content: {round(final_gc, 4)}", normal_style))
elements.append(Paragraph(f"GC3 content: {round(final_gc3, 4)}", normal_style))
elements.append(Spacer(1, 0.2 * inch))

elements.append(Paragraph("2.1 Final Optimized Sequence", sub_heading_style))
fasta_text = f">{record.id}\n{format_fasta(final_seq)}"
elements.append(Preformatted(fasta_text, fasta_style))

elements.append(PageBreak())
elements.append(Paragraph("2.2 Codon Usage Comparison", sub_heading_style))
grid_images = sorted([f for f in os.listdir(grid_dir) if f.endswith(".png")])
for img in grid_images:
    img_path = os.path.join(grid_dir, img)
    elements.append(Image(img_path, width=6.5 * inch, height=9 * inch))
    elements.append(PageBreak())

#############################################
# Existing Wright plot
#############################################

elements.append(Paragraph("3. Genome Codon Bias (Wright Plot)", heading_style))
elements.append(Spacer(1, 0.1 * inch))
elements.append(Image(wright_plot_file, width=5 * inch, height=4 * inch))

#############################################
# New diagnostics section
#############################################

elements.append(PageBreak())
elements.append(Paragraph("4. Extended Diagnostic Plots", heading_style))
elements.append(Spacer(1, 0.1 * inch))

plot_notes = [
    (
        "enc_distribution.png",
        "4.1 ENC Distribution",
        "Shows effective number of codons (ENC) across CDS. Lower ENC (<35) suggests stronger codon bias; values >50 indicate weaker bias.",
    ),
    (
        "gc3_distribution.png",
        "4.2 GC3 Distribution",
        "Distribution of third-position GC. Broad spread suggests heterogeneous mutational pressure; many genes above 0.65 indicate GC-rich codon endings.",
    ),
    (
        "wright_hexbin.png",
        "4.3 Wright Density Plot",
        "ENC vs GC3 density relative to Wright expected curve. Genes far below expected line indicate selection beyond mutational bias.",
    ),
    (
        "wright_residuals.png",
        "4.4 Wright Residuals",
        "Residuals (ENC observed - expected). Strongly negative residuals (e.g., < -5) suggest highly selected codon usage.",
    ),
    (
        "rscu_heatmap.png",
        "4.5 RSCU Heatmap",
        "Relative synonymous codon usage (RSCU) for genome, ribosomal, and final sequence. RSCU >1.5 indicates codon enrichment in a context.",
    ),
    (
        "cai_distributions.png",
        "4.6 CAI Distributions",
        "Genome CDS CAI distributions under genome/ribosomal weights, with final sequence marked. Final CAI above upper quartile implies strong adaptation.",
    ),
    (
        "sliding_gc_profile.png",
        "4.7 Sliding-window GC",
        "Local GC profile against target and ±1σ bands. Frequent excursions beyond ±1σ (0.05) highlight potential local synthesis/expression risk.",
    ),
    (
        "rare_codon_runs.png",
        "4.8 Rare-codon Run Lengths",
        "Histogram of consecutive rare codon runs (here rare ~ weight <0.1). Runs >=3 can be practical warning zones for translation slowdown.",
    ),
    (
        "pareto_scatter_matrix.png",
        "4.9 Pareto Objective Scatter Matrix",
        "Pairwise trade-offs among CAI_r, CAI_g, GC score, and splice score. A wider front indicates richer optimization choices.",
    ),
    (
        "candidate_diversity.png",
        "4.10 Candidate Diversity",
        "Pairwise codon Hamming distance among exported candidates. Very low diversity (<0.05) may imply over-convergence.",
    ),
    (
        "ga_objective_trajectories.png",
        "4.11 GA Objective/Control Trajectories",
        "Tracks fitness dynamics and controls (e.g., mutation rate, front size). Flat best-fitness with shrinking front can indicate stagnation.",
    ),
    (
        "synonymous_entropy.png",
        "4.12 Synonymous Entropy",
        "Per-amino-acid normalized synonymous entropy (0..1). Values near 0 imply codon collapse; near 1 imply balanced synonymous usage.",
    ),
    (
        "codon_usage_pca.png",
        "4.13 Codon Usage PCA",
        "Projects CDS codon-usage vectors with final sequence overlay. Final point far from CDS cloud suggests less genome-like codon composition.",
    ),
    (
        "motif_risk_map.png",
        "4.14 Motif Risk Map",
        "Positions of polyA motifs, donor motifs, and GT...AG spans. Dense motif clusters are practical warning regions for transcript safety.",
    ),
]

for fname, title, desc in plot_notes:
    fpath = os.path.join(diagnostics_dir, fname)
    if os.path.exists(fpath):
        elements.append(Paragraph(title, sub_heading_style))
        elements.append(Paragraph(desc, normal_style))
        elements.append(Spacer(1, 0.08 * inch))
        elements.append(Image(fpath, width=6.4 * inch, height=3.2 * inch))
        elements.append(Spacer(1, 0.18 * inch))

#############################################
# Interpretation Notes
#############################################

elements.append(PageBreak())
elements.append(Paragraph("5. Interpretation Notes", heading_style))
elements.append(Spacer(1, 0.2 * inch))

notes = [
    "Optimization balances translational efficiency and genome-wide codon conformity.",
    "GC content is constrained using a Gaussian penalty centered on genome average.",
    "Cryptic splice-site penalties were applied during genetic algorithm optimization.",
    "Use diagnostics jointly: one favorable metric alone is not sufficient for sequence selection.",
    "Final sequence should be validated experimentally before synthesis.",
]

elements.append(ListFlowable([ListItem(Paragraph(n, normal_style)) for n in notes], bulletType="bullet"))

doc.build(elements)
