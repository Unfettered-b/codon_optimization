#!/usr/bin/env python3

import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

genome_fasta = snakemake.input.genome
gff_file = snakemake.input.gff
output_fasta = snakemake.output[0]

# Load genome sequences
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Store CDS features grouped by gene (Parent attribute)
cds_dict = defaultdict(list)

with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) != 9:
            continue
        
        seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
        
        if feature_type != "CDS":
            continue
        
        start = int(start)
        end = int(end)
        
        # Parse attributes
        attr_dict = {}
        for item in attributes.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attr_dict[key] = value
        
        parent_id = attr_dict.get("Parent")
        if not parent_id:
            continue
        
        cds_dict[parent_id].append({
            "seqid": seqid,
            "start": start,
            "end": end,
            "strand": strand
        })

# Extract and write CDS sequences
with open(output_fasta, "w") as out:
    for gene_id, features in cds_dict.items():
        
        # Sort exons by genomic position
        features_sorted = sorted(features, key=lambda x: x["start"])
        
        seq_fragments = []
        strand = features_sorted[0]["strand"]
        
        for feat in features_sorted:
            chrom_seq = genome[feat["seqid"]].seq
            fragment = chrom_seq[feat["start"] - 1 : feat["end"]]
            seq_fragments.append(fragment)
        
        cds_seq = Seq("").join(seq_fragments)
        
        # Reverse complement if on negative strand
        if strand == "-":
            cds_seq = cds_seq.reverse_complement()
        
        # Skip sequences not divisible by 3
        if len(cds_seq) % 3 != 0:
            continue
        
        out.write(f">{gene_id}\n{cds_seq}\n")