#!/usr/bin/env python3

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: python replace_leucine_with_ctc.py input.fasta output.fasta")
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

# All leucine codons
leu_codons = {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}

with open(output_fasta, "w") as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()

        # Ensure sequence length is multiple of 3
        if len(seq) % 3 != 0:
            raise ValueError(f"Sequence {record.id} is not divisible by 3")

        new_seq = []

        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]

            if codon in leu_codons:
                new_seq.append("CTC")
            else:
                new_seq.append(codon)

        modified_sequence = "".join(new_seq)

        out_handle.write(f">{record.id}\n{modified_sequence}\n")