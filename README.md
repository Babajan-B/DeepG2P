# DeepG2P-Celiac

Code for ranking celiac disease candidate genes from public GWAS Catalog records and testing the ranking
in independent summary statistics (IEU OpenGWAS ieu-a-276 and FinnGen release 12).

This repository accompanies the manuscript *"DeepG2P-Celiac: multi-source evidence synthesis for ranking
celiac disease risk genes"* (Shaik NA, Elango R, Awan Z, Banaganapalli B; submitted). The scripts in
`scripts/` regenerate every number, table and figure in that manuscript from the downloaded public files.

No data are stored in this repository. All inputs are public summary-level data that you download yourself
(see below), and no individual-level genotype data were used for the manuscript.

## What the analysis does

1. Extracts celiac disease associations from the NHGRI-EBI GWAS Catalog (259 associations, 16 studies,
   release indexed 18 May 2026).
2. Splits each association into one record per catalog-mapped gene and scores each gene as
   `ES = min(-log10 p_best, 300) + n_records`. The extended MHC (chr6:25–34 Mb, GRCh38) is reported
   separately.
3. Compares the score with ranking by p value alone, and runs a sensitivity analysis with a seven-edge
   curated gene graph.
4. Looks up a 10-SNP candidate panel (`configs/celiac_snp_panel.csv`) in ieu-a-276 and FinnGen R12.
5. Tests whether the score predicts gene-level replication in FinnGen (AUC, Spearman correlation,
   cluster bootstrap) and checks risk-allele direction.

## Installation

Requires Python 3.10 or later. The FinnGen step also needs `tabix` from
[HTSlib](https://www.htslib.org/) built with network (libcurl) support.

```bash
git clone https://github.com/Babajan-B/DeepG2P.git
cd DeepG2P
python -m venv .venv
source .venv/bin/activate
pip install -e ".[analysis,dev]"
pytest
```

## Reproducing the manuscript

Run every command from the repository root.

**1. Download the GWAS Catalog association file** (about 60 MB):

```bash
mkdir -p data/raw/open_source
curl -L -o data/raw/open_source/gwas-catalog-associations-full.zip \
  https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip
python scripts/extract_gwas_catalog.py
```

The catalog is updated continuously, so a later release may return more than the 259 associations used
in the manuscript.

**2. Download IEU OpenGWAS dataset ieu-a-276** (Dubois et al. 2010; GWAS-VCF, GRCh37) from
<https://gwas.mrcieu.ac.uk/datasets/ieu-a-276/> (a free OpenGWAS account is required) and save it as
`data/raw/open_source/ieu-a-276/ieu-a-276.vcf.gz`. Then extract the top variants:

```bash
deepg2p-ieu-vcf \
  --vcf data/raw/open_source/ieu-a-276/ieu-a-276.vcf.gz \
  --output data/processed/open_source_evidence/ieu-a-276/ieu_a_276_top1000_summary_stats.csv
```

**3. Gene ranking, panel look-up, figures 1–4 and supplementary tables:**

```bash
python scripts/make_figures.py
```

**4. FinnGen R12 replication and figure 5.** FinnGen summary statistics for `K11_COELIAC` are queried
remotely with tabix from the public FinnGen bucket, so nothing large is downloaded. Panel SNP positions on
GRCh38 are fetched from the Ensembl REST API.

```bash
python scripts/replicate_finngen.py
```

Outputs are written to `results/`:

| Folder | Contents |
|---|---|
| `results/figures/` | Figures 1–5 (PNG and TIFF) |
| `results/supplementary/` | Gene scores (non-MHC and MHC), top ieu-a-276 variants, supplementary figure |
| `results/audit/` | Derived numbers, panel look-ups, FinnGen replication tables and direction check |

## Repository layout

```
configs/    candidate SNP panel and pipeline settings
docs/       data policy and workflow notes
scripts/    manuscript analysis scripts (extraction, scoring, figures, FinnGen replication)
src/        the deepg2p_celiac Python package
tests/      unit tests (pytest)
data/       empty placeholders; downloaded data go here and are ignored by git
results/    empty placeholder; outputs go here and are ignored by git
```

## Other modules

The `deepg2p_celiac` package also contains modules for individual-level data: SNP-panel extraction from
consumer genotype files, phenotype cleaning, PLINK QC command building, a baseline classifier, a
graph-smoothed gene model and a small FastAPI service. These modules need real, consented genotype and
phenotype files (see [docs/data_policy.md](docs/data_policy.md)). **They were not used in the manuscript,**
and no trained model contributed to its results. Their command-line entry points are listed in
`pyproject.toml`, and [docs/final_workflow.md](docs/final_workflow.md) describes their use.

## Limitations

- Gene assignment uses the GWAS Catalog `MAPPED_GENE` field. The ranking does not identify causal genes.
- Record counts track how often a locus was published, and SNPs are not clumped by linkage disequilibrium.
- Nearly all of the input data come from people of European ancestry.
- HLA tag SNPs were designed for haplotype calling and cannot be evaluated from single-SNP summary
  statistics.

This software is for research only and is not intended for clinical diagnosis.

## Data sources and acknowledgements

- NHGRI-EBI GWAS Catalog: Sollis E et al. (2023) *Nucleic Acids Res* 51:D977–D985.
- IEU OpenGWAS: Elsworth B et al. (2020) *bioRxiv*. Dataset ieu-a-276: Dubois PCA et al. (2010)
  *Nat Genet* 42:295–302.
- FinnGen: Kurki MI et al. (2023) *Nature* 613:508–518.

Each data source has its own terms of use, and they apply when you download the data.

This work was funded by the Deanship of Scientific Research (DSR) at King Abdulaziz University, Jeddah,
Saudi Arabia, under grant No. (G: 444-140-1442).

## Authors

- Noor Ahmad Shaik, Department of Genetic Medicine, King Abdulaziz University
  ([ORCID 0000-0001-6674-5409](https://orcid.org/0000-0001-6674-5409))
- Ramu Elango, Department of Genetic Medicine, King Abdulaziz University
  ([ORCID 0000-0003-4264-2135](https://orcid.org/0000-0003-4264-2135))
- Zuhier Awan, Department of Clinical Biochemistry, King Abdulaziz University
  ([ORCID 0000-0003-2850-9047](https://orcid.org/0000-0003-2850-9047))
- Babajan Banaganapalli (corresponding author, bbabajan@kau.edu.sa), Department of Genetic Medicine,
  King Abdulaziz University ([ORCID 0000-0001-8268-3412](https://orcid.org/0000-0001-8268-3412))

## Citation

If you use this code, please cite the repository (see [CITATION.cff](CITATION.cff)) and, once it is
published, the accompanying article.

## License

Released under the MIT License. See [LICENSE](LICENSE).
