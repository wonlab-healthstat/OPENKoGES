# OpenKoGES Analysis Code

This repository contains all analysis code used in the following manuscript:

**The Korean Genome and Epidemiology Study: A Large-Scale Genomic and Longitudinal Resource for Population-Scale Analysis**

## Authors of the Code

The code in this repository was written by **sjim** and **joh**.

For any inquiries, please contact **soojinim1@gmail.com**.

## Directory Overview

### `CODE_manuscript_final/`

Analysis scripts used to generate the manuscript figures and tables:

| Path | Description |
|------|-------------|
| `1_Expansion_of_the_KoGES_resource/` | KoGES resource expansion (Figure 1, Extended Data 1, Supple Figure/Table 1) |
| `2_Population_structure_relatedness_allele_frequency/` | Population structure, relatedness, and allele frequency analyses (Figure 2) |
| `3_GWAS/` | GWAS analyses (BOLT-LMM / SAIGE), lead variant annotation, heritability, and genetic correlation (Figure 3) |
| `4_70K_Comparison/` | Comparison with the ~70K cohort (Supple Figure 6, Supple Table 5) |
| `5_Novelty_Replication/` | Novelty assessment, replication in EAS/EUR cohorts, COJO, and winner’s curse analyses |
| `6_Heritability_GC/` | Heritability and genetic correlation summaries (Figure 5) |
| `ExtendedData_export/` | Scripts to export Extended Data tables for publication |

### `RESULT_manuscript_final/`

Output files corresponding to manuscript figures and tables:

- `Figure1/`–`Figure5/` — main figures
- `ExtendedData/` — Extended Data outputs
- `SuppleFigure*/` — supplementary figures
- `SuppleTable*/` — supplementary tables
