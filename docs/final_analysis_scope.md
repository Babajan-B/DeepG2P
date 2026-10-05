# Final Analysis Scope

DeepG2P-Celiac can be completed in two scientifically honest modes.

## Mode 1: Individual-Level Validation

Requires real celiac case/control genotype and phenotype files.

Outputs:

- Target SNP matrix.
- QC matrix.
- Baseline model validation.
- Graph-smoothed model validation.
- Cross-validation metrics.
- Calibration.
- Permutation testing.
- Feature attribution.

## Mode 2: Open-Source Evidence Synthesis

Does not require individual-level genotype files.

Inputs:

- Public GWAS Catalog-style association tables.
- Open Targets-style variant/gene/trait evidence tables.
- Any curated celiac association CSV with variant, gene, trait, and p-value fields.

Outputs:

- Normalized celiac association table.
- Gene evidence scores.
- Graph-context gene ranking.
- JSON evidence report.

## Boundary

Mode 2 is not predictive validation. It is a literature and open-association evidence layer. It is suitable for ending the project as an analysis framework when real individual-level celiac genotype validation data is not available.

