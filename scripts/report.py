#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image,
    ListFlowable, ListItem
)
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet

#############################################
# Inputs
#############################################

ga_log_file = snakemake.input.ga_log
final_seq_file = snakemake.input.final_seq
wright_plot_file = snakemake.input.wright

output_pdf = snakemake.output[0]

#############################################
# Load Data
#############################################

ga_log_df = pd.read_csv(ga_log_file, sep="\t")

record = next(SeqIO.parse(final_seq_file, "fasta"))
final_seq = str(record.seq)

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
# Final Sequence Metrics
#############################################

elements.append(Paragraph("2. Final Optimized Sequence Metrics", heading_style))
elements.append(Spacer(1, 0.2 * inch))

elements.append(Paragraph(f"Sequence length: {len(final_seq)} nt", normal_style))
elements.append(Paragraph(f"GC content: {round(final_gc, 4)}", normal_style))
elements.append(Paragraph(f"GC3 content: {round(final_gc3, 4)}", normal_style))
elements.append(Spacer(1, 0.4 * inch))

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