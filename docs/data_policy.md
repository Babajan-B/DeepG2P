# Data Policy

This project must use real, consented, public, or user-provided genotype and phenotype files.

## Not Allowed

- Synthetic genotype datasets for model training.
- Randomly generated case/control labels.
- Simulated phenotype files presented as real data.
- Full public genotype archive downloads on this machine.
- Clinical diagnosis claims.

## Allowed

- Real consumer genotype text files with consent or valid public licensing.
- Public phenotype metadata with clear provenance.
- Small in-memory test fixtures used only for software tests.
- Reduced SNP-panel matrices extracted from real raw genotype files.

## Required Provenance

Every real data run should document:

- Source name.
- Download or access date.
- License or consent basis.
- Number of raw genotype files processed.
- Number of samples after phenotype cleaning.
- Number of SNPs extracted.

