#!/usr/bin/env python3

import re

gff_file = snakemake.input.gff
output_file = snakemake.output[0]

# Keywords indicating ribosomal / translation genes
ribosomal_patterns = [
    r"ribosomal protein",
    r"rpl\d+",
    r"rps\d+",
    r"60s ribosomal",
    r"40s ribosomal",
    r"elongation factor",
    r"translation initiation factor",
    r"translation elongation factor"
]

compiled_patterns = [re.compile(p, re.IGNORECASE) for p in ribosomal_patterns]

ribosomal_ids = set()

with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) != 9:
            continue

        seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

        # We check mRNA or CDS entries
        if feature_type not in ["mRNA", "transcript", "CDS"]:
            continue

        attr_dict = {}
        for item in attributes.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attr_dict[key] = value

        # Search in product / description fields
        annotation_fields = [
            attr_dict.get("product", ""),
            attr_dict.get("Name", ""),
            attr_dict.get("gene", ""),
            attr_dict.get("description", "")
        ]

        annotation_text = " ".join(annotation_fields)

        for pattern in compiled_patterns:
            if pattern.search(annotation_text):
                # Prefer transcript ID if available
                gene_id = attr_dict.get("ID") or attr_dict.get("Parent")
                if gene_id:
                    ribosomal_ids.add(gene_id)
                break

# Write unique IDs
with open(output_file, "w") as out:
    for rid in sorted(ribosomal_ids):
        out.write(f"{rid}\n")