#!/usr/bin/env python3

import pandas as pd
import re
import os
from Bio import SeqIO
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image,
    ListFlowable, ListItem, Preformatted, PageBreak
)
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors


#############################################
# Inputs
#############################################

ga_log_file = snakemake.input.ga_log
final_seq_file = snakemake.input.final_seq
wright_plot_file = snakemake.input.wright
fitness_plot_file = snakemake.input.fitness_plot
cai_dir = snakemake.input.cai_dir
grid_dir = snakemake.input.codon_usage_grid

output_pdf = snakemake.output[0]

#############################################
# Load Data
#############################################

ga_log_df = pd.read_csv(ga_log_file, sep="\t")

record = next(SeqIO.parse(final_seq_file, "fasta"))
final_seq = str(record.seq)


fasta_style = ParagraphStyle(
    name="FastaStyle",
    fontName="Courier",
    fontSize=8,
    leading=10
)

def format_fasta(seq):
    lines = [seq[i:i+60] for i in range(0, len(seq), 60)]
    return "\n".join(lines)

def extract_generation(filename):
    match = re.search(r"generation_(\d+)", filename)
    if match:
        return int(match.group(1))
    return -1

#############################################
# Compute Final Metrics
#############################################

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
elements.append(Spacer(1, 0.4 * inch))

#############################################
# GA Progression plot
#############################################

elements.append(Paragraph("1.1 Fitness Progression", sub_heading_style))
elements.append(Spacer(1, 0.2 * inch))
elements.append(Image(fitness_plot_file, width=5*inch, height=3*inch))
elements.append(Spacer(1, 0.4 * inch))

#############################################
# CAI Snapshot Plots
#############################################

elements.append(PageBreak())
elements.append(Paragraph("1.2 CAI Positional Plots Across Generations", styles["Heading2"]))
elements.append(Spacer(1, 0.3 * inch))

# Get and sort images
cai_images = [
    f for f in os.listdir(cai_dir)
    if f.endswith(".png")
]

cai_images.sort(key=extract_generation)

for img_file in cai_images:

    gen_number = extract_generation(img_file)
    img_path = os.path.join(cai_dir, img_file)

    elements.append(Paragraph(f"Generation {gen_number}", styles["Heading3"]))
    elements.append(Spacer(1, 0.1 * inch))

    elements.append(Image(img_path, width=6*inch, height=1.5*inch))
    elements.append(Spacer(1, 0.3 * inch))

#############################################
# Final Sequence Metrics
#############################################

elements.append(Paragraph("2. Final Optimized Sequence Metrics", heading_style))
elements.append(Spacer(1, 0.2 * inch))

elements.append(Paragraph(f"Sequence length: {len(final_seq)} nt", normal_style))
elements.append(Paragraph(f"GC content: {round(final_gc, 4)}", normal_style))
elements.append(Paragraph(f"GC3 content: {round(final_gc3, 4)}", normal_style))
elements.append(Spacer(1, 0.4 * inch))


elements.append(Paragraph("2.1 Final Optimized Sequence", sub_heading_style))
elements.append(Spacer(1, 0.2 * inch))

fasta_text = f">{record.id}\n{format_fasta(final_seq)}"
elements.append(Preformatted(fasta_text, fasta_style))

elements.append(PageBreak())
elements.append(Paragraph("2.2 Codon Usage Comparison", sub_heading_style))

grid_images = sorted(
    [f for f in os.listdir(grid_dir) if f.endswith(".png")]
)

for img in grid_images:
    img_path = os.path.join(grid_dir, img)
    elements.append(Image(img_path, width=6.5*inch, height=9*inch))
    elements.append(PageBreak())

#############################################
# Wright Plot
#############################################

elements.append(Paragraph("3. Genome Codon Bias (Wright Plot)", heading_style))
elements.append(Spacer(1, 0.2 * inch))

img = Image(wright_plot_file, width=5*inch, height=4*inch)
elements.append(img)
elements.append(Spacer(1, 0.4 * inch))

#############################################
# Interpretation Notes
#############################################

elements.append(Paragraph("4. Interpretation Notes", heading_style))
elements.append(Spacer(1, 0.2 * inch))

notes = [
    "Optimization balances translational efficiency and genome-wide codon conformity.",
    "GC content is constrained using a Gaussian penalty centered on genome average.",
    "Cryptic splice-site penalties were applied during genetic algorithm optimization.",
    "Final sequence should be validated experimentally before synthesis."
]

elements.append(ListFlowable(
    [ListItem(Paragraph(n, normal_style)) for n in notes],
    bulletType='bullet'
))

#############################################
# Build PDF
#############################################

doc.build(elements)