#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import PageBreak
from reportlab.platypus import ListFlowable, ListItem

#############################################
# Inputs
#############################################

native_metrics_file = snakemake.input.native
ga_log_file = snakemake.input.ga_log
final_seq_file = snakemake.input.final_seq
wright_plot_file = snakemake.input.wright

output_pdf = snakemake.output[0]

#############################################
# Load Data
#############################################

native_df = pd.read_csv(native_metrics_file, sep="\t")
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

#############################################
# Title
#############################################

elements.append(Paragraph("Codon Optimization Report", title_style))
elements.append(Spacer(1, 0.3 * inch))

#############################################
# Native Metrics Section
#############################################

elements.append(Paragraph("1. Baseline (Native) Metrics", styles["Heading2"]))
elements.append(Spacer(1, 0.2 * inch))

native_data = [["Gene ID", "Length (nt)", "CAI_ribo", "CAI_genome", "Composite", "GC", "GC3"]]

for _, row in native_df.iterrows():
    native_data.append([
        row.get("gene_id", ""),
        row.get("length_nt", ""),
        row.get("CAI_ribo", ""),
        row.get("CAI_genome", ""),
        row.get("Composite_score", ""),
        row.get("GC", ""),
        row.get("GC3", "")
    ])

table = Table(native_data, repeatRows=1)
table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
]))
elements.append(table)
elements.append(Spacer(1, 0.4 * inch))

#############################################
# GA Summary
#############################################

elements.append(Paragraph("2. Genetic Algorithm Optimization Summary", styles["Heading2"]))
elements.append(Spacer(1, 0.2 * inch))

elements.append(Paragraph(f"Total generations: {final_generation}", normal_style))
elements.append(Paragraph(f"Best composite fitness achieved: {round(best_fitness, 5)}", normal_style))
elements.append(Spacer(1, 0.4 * inch))

#############################################
# Final Sequence Metrics
#############################################

elements.append(Paragraph("3. Final Optimized Sequence Metrics", styles["Heading2"]))
elements.append(Spacer(1, 0.2 * inch))

elements.append(Paragraph(f"Sequence length: {len(final_seq)} nt", normal_style))
elements.append(Paragraph(f"GC content: {round(final_gc, 4)}", normal_style))
elements.append(Paragraph(f"GC3 content: {round(final_gc3, 4)}", normal_style))
elements.append(Spacer(1, 0.4 * inch))

#############################################
# Wright Plot
#############################################

elements.append(Paragraph("4. Genome Codon Bias (Wright Plot)", styles["Heading2"]))
elements.append(Spacer(1, 0.2 * inch))

img = Image(wright_plot_file, width=5*inch, height=4*inch)
elements.append(img)
elements.append(Spacer(1, 0.4 * inch))

#############################################
# Interpretation Notes
#############################################

elements.append(Paragraph("5. Interpretation Notes", styles["Heading2"]))
elements.append(Spacer(1, 0.2 * inch))

notes = [
    "CAI_ribo reflects adaptation to highly expressed ribosomal genes.",
    "CAI_genome reflects conformity to genome-wide codon bias.",
    "Composite fitness balances translational efficiency and genomic realism.",
    "GC content is controlled to remain near genome average.",
    "Cryptic splice site penalties were applied during optimization."
]

elements.append(ListFlowable(
    [ListItem(Paragraph(n, normal_style)) for n in notes],
    bulletType='bullet'
))

#############################################
# Build PDF
#############################################

doc.build(elements)