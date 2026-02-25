# 🧬 Fungal Codon Optimization Pipeline

### Multi-Objective CAI-Driven Genetic Algorithm for *Aureobasidium pullulans*

------------------------------------------------------------------------

![Snakemake](https://img.shields.io/badge/workflow-Snakemake-blue)
![Python](https://img.shields.io/badge/python-3.10%2B-green)
![R](https://img.shields.io/badge/R-coRdon-orange)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

A reproducible **Snakemake-based codon optimization framework** designed
for fungal nuclear expression systems.

This pipeline integrates:

-   📊 Genome-wide codon bias analysis\
-   🧬 Dual CAI model construction (ribosomal + genome)\
-   🧠 Multi-objective genetic algorithm optimization\
-   📈 GC balancing with Gaussian constraint\
-   ⚠️ Cryptic splice-site suppression\
-   📄 Automated PDF reporting

Designed specifically for ***Aureobasidium pullulans***, but adaptable
to other fungi.

------------------------------------------------------------------------

# ✨ Why This Pipeline?

Codon optimization is often oversimplified to:

> "Maximize CAI."

This pipeline instead balances:

-   Translational efficiency\
-   Genome realism\
-   Mutation bias\
-   Transcript safety

Optimization objective:

F = (CAI_ribo\^α) · (CAI_genome\^β) · P_GC · P_splice

This avoids:

-   GC over-inflation\
-   Artificial codon collapse\
-   Cryptic intron formation\
-   Premature polyadenylation

------------------------------------------------------------------------

# 📂 Repository Structure

    .
    ├── Snakefile
    ├── config.yaml
    ├── data/
    ├── scripts/
    ├── results/
    └── README.md

------------------------------------------------------------------------

# 🔬 Pipeline Workflow

## 1️⃣ Genome Codon Bias Analysis

-   Extract CDS\
-   Compute ENC and GC3\
-   Generate Wright plot\
-   Identify ribosomal genes

Outputs: - `all_cds.fasta` - `gene_metrics.tsv` - `wright_plot.png` -
`ribosomal_ids.txt`

------------------------------------------------------------------------

## 2️⃣ CAI Model Construction

Two reference models:

  Model           Source               Purpose
  --------------- -------------------- -----------------------------
  Ribosomal CAI   Ribosomal CDS only   High-expression reference
  Genome CAI      All CDS              Genome conformity reference

Outputs: - `codon_weights_ribo.json` - `codon_weights_genome.json`

------------------------------------------------------------------------

## 3️⃣ Genetic Algorithm Optimization

Input: - Protein sequence (FASTA)

Optimization components:

✔ CAI adaptation to ribosomal genes\
✔ Genome-wide conformity\
✔ Gaussian GC constraint\
✔ Cryptic splice suppression

Outputs: - `optimized_sequence.fasta` - `ga_log.tsv`

------------------------------------------------------------------------

## 4️⃣ Fungal Transcript Safety Filter

Final validation checks:

-   Internal stop codons\
-   PolyA signals (AATAAA, ATTAAA)\
-   Strong donor motifs\
-   GT...AG intron-like patterns\
-   Extreme GC windows

Outputs: - `final_sequence.fasta` - `fungal_warnings.txt`

------------------------------------------------------------------------

## 5️⃣ Automated Report

Generates:

-   Baseline metrics\
-   GA convergence summary\
-   Final GC statistics\
-   Wright plot\
-   Interpretation notes

Output: - `optimization_report.pdf`

------------------------------------------------------------------------

# ⚙️ Configuration

All parameters are controlled via `config.yaml`.

Example:

``` yaml
ga:
  population_size: 120
  generations: 250
  mutation_rate: 0.05
  elite_size: 5

fitness:
  alpha: 0.7
  beta: 0.3

gc:
  target: 0.52
  sigma: 0.05

snapshot:
  count: 10

fungal:
  gc_window: 50
  low_gc: 0.30
  high_gc: 0.70
```

------------------------------------------------------------------------

# 🚀 Running the Pipeline

Run full workflow:

``` bash
snakemake --use-conda --cores 8
```

Run specific target:

``` bash
snakemake results/optimized_sequence.fasta
```

Generate DAG:

``` bash
snakemake --dag | dot -Tpdf > dag.pdf
```

------------------------------------------------------------------------

# 🧪 Requirements

-   Python ≥ 3.10\
-   Snakemake\
-   Biopython\
-   pandas\
-   R + coRdon\
-   reportlab

Conda environments recommended.

------------------------------------------------------------------------

# 🔬 Applications

-   Recombinant protein expression in fungi\
-   Synthetic biology\
-   Codon engineering research\
-   Expression vector design

------------------------------------------------------------------------

# 📌 Notes

-   Assumes nuclear expression.\
-   Kozak context handled at vector level.\
-   RNA structure penalties can be added in future versions.

------------------------------------------------------------------------

# 👨‍🔬 Author

Developed for fungal codon engineering and translational optimization
research.
