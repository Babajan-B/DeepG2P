# Final Workflow

This project is complete as an analysis-first celiac genetics framework.

## Preferred Completion Mode: Open-Source Evidence

Use this mode when individual-level genotype and phenotype data are not available.

1. Download a public GWAS Catalog/Open Targets-style celiac association CSV.
2. Run:

```bash
deepg2p-open-evidence \
  --associations path/to/association_file.csv \
  --output-dir data/processed/open_source_evidence \
  --trait-filter celiac
```

3. Create candidate panel updates:

```bash
deepg2p-update-panel \
  --panel configs/celiac_snp_panel.csv \
  --evidence-scores data/processed/open_source_evidence/gene_evidence_scores.csv \
  --output data/processed/open_source_evidence/celiac_snp_panel_candidates.csv
```

4. Render the final report:

```bash
deepg2p-final-report \
  --evidence-report data/processed/open_source_evidence/open_source_evidence_report.json \
  --evidence-scores data/processed/open_source_evidence/gene_evidence_scores.csv \
  --output data/processed/open_source_evidence/final_analysis_report.html
```

## Optional Mode: Individual-Level Validation

Use this mode only when real case/control genotype files are available.

```bash
deepg2p-run-pipeline \
  --genotype-dir data/raw/genotypes \
  --phenotypes data/raw/phenotypes.csv \
  --panel configs/celiac_snp_panel.csv \
  --output-dir data/processed/run \
  --max-gb 6
```

## Boundary

Open-source evidence mode is not predictive validation. It is a gene/locus prioritization analysis based on published association evidence.

