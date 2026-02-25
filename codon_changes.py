#!/usr/bin/env python3

from Bio.Seq import Seq

from collections import Counter

def analyze_codon_changes(seq1, seq2):
    """
    Returns:
    (
        total_codons,
        identical_codons,
        changed_codons,
        synonymous_changes,
        nonsynonymous_changes,
        codon_change_frequency_dict
    )
    """

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length.")

    if len(seq1) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3.")

    total_codons = len(seq1) // 3
    identical_codons = 0
    changed_codons = 0
    synonymous = 0
    nonsynonymous = 0

    codon_changes = []

    for i in range(0, len(seq1), 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]

        if codon1 == codon2:
            identical_codons += 1
        else:
            changed_codons += 1
            codon_changes.append((codon1, codon2))

            aa1 = str(Seq(codon1).translate())
            aa2 = str(Seq(codon2).translate())

            if aa1 == aa2:
                synonymous += 1
            else:
                nonsynonymous += 1

    codon_freq = Counter(codon_changes)

    return (
        total_codons,
        identical_codons,
        changed_codons,
        synonymous,
        nonsynonymous,
        codon_freq
    )

# ================= Example usage =================

if __name__ == "__main__":
    seq1 = """ATGGTCTCTAAGGGAGAGGCAGTCATCAAGGAATTCATGCGCTTCAAAGTTCATATGGAGGGCTCAATGAACGGCCACGAATTCGAAATTGAAGGAGAAGGCGAGGGCAGACCTTACGAGGGCACACAAACTGCAAAACTAAAAGTTACTAAAGGAGGCCCTCTCCCTTTTTCATGGGATATCCTTTCACCTCAATTCATGTACGGCAGCAGAGCTTTCATCAAGCATCCTGCCGACATCCCCGATTACTATAAGCAATCTTTCCCTGAAGGCTTCAAATGGGAGAGAGTCATGAACTTCGAAGATGGCGGAGCTGTGACCGTCACTCAAGATACCAGCCTTGAGGATGGCACACTTATCTATAAAGTCAAGCTCCGCGGCACGAACTTTCCACCTGATGGCCCTGTTATGCAAAAGAAGACTATGGGCTGGGAAGCTAGCACCGAGCGCCTTTACCCTGAGGATGGAGTCCTTAAGGGCGACATTAAAATGGCACTCCGCCTCAAGGACGGCGGCAGATACTTGGCAGATTTCAAAACAACATATAAGGCTAAGAAGCCTGTTCAAATGCCTGGAGCTTACAACGTTGATAGAAAGCTTGACATCACCTCGCATAACGAGGATTACACTGTCGTCGAGCAATACGAGAGATCAGAGGGCCGACACAGCACCGGAGGAATGGATGAACTCTACAAGCTTGAAGGCTCG""".replace("\n", "") # replace with full seq
    seq2 = """ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCTCTCTCCCAGCGACACATGAGTTACACATCTTTGGCTCCATCAACGGTGTGGACTTTGACATGGTGGGTCAGGGCACCGGCAATCCAAATGATGGTTATGAGGAGTTAAACCTGAAGTCCACCAAGGGTGACCTCCAGTTCTCCCCCTGGATTCTGGTCCCTCATATCGGGTATGGCTTCCATCAGTACCTGCCCTACCCTGACGGGATGTCGCCTTTCCAGGCCGCCATGGTAGATGGCTCCGGCTACCAAGTCCATCGCACAATGCAGTTTGAAGATGGTGCCTCCCTTACTGTTAACTACCGCTACACCTACGAGGGAAGCCACATCAAAGGAGAGGCCCAGGTGAAGGGGACTGGTTTCCCTGCTGACGGTCCTGTGATGACCAACTCGCTGACCGCTGCGGACTGGTGCAGGTCGAAGAAGACTTACCCCAACGACAAAACCATCATCAGTACCTTTAAGTGGAGTTACACCACTGGAAATGGCAAGCGCTACCGGAGCACTGCGCGGACCACCTACACCTTTGCCAAGCCAATGGCGGCTAACTATCTGAAGAACCAGCCGATGTACGTGTTCCGTAAGACGGAGCTCAAGCACTCCAAGACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATGTGATGGGCATGGACGAGCTGTACAAG""".replace("\n", "")  # replace with full seq

    result = analyze_codon_changes(seq1, seq2)
    print("Results:", result)