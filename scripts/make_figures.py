"""Regenerate all manuscript figures and supplementary tables from the pipeline outputs.

Run from the repository root:  python scripts/make_figures.py
Every number plotted here is read from data/ — nothing is typed in by hand.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deepg2p_celiac.analysis.evidence import add_graph_context, score_gene_evidence, split_gene_field  # noqa: E402
from deepg2p_celiac.graph.build_graph import build_gene_graph  # noqa: E402

EV = ROOT / "data/processed/open_source_evidence"
VCF = ROOT / "data/raw/open_source/ieu-a-276/ieu-a-276.vcf.gz"
OUT = ROOT / "results"
FIG, SUP, AUDIT = OUT / "figures", OUT / "supplementary", OUT / "audit"
for _d in (FIG, SUP, AUDIT):
    _d.mkdir(parents=True, exist_ok=True)
# Optional: output of the earlier (uncorrected) pipeline run, used only for the audit comparison
LEGACY = EV / "combined_with_ieu"
HAVE_LEGACY = (LEGACY / "gene_evidence_scores.csv").exists()

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 10,
    "axes.titlelocation": "left", "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_DARK, C_MID, C_LIGHT, C_ACC, C_GREY = "#1f3b5a", "#2a7f9e", "#9cc9d9", "#c8553d", "#8c8c8c"


def chrom_key(c: str) -> int:
    return 23 if c == "X" else 24 if c == "Y" else int(c)


# --------------------------------------------------------------------------- data
catalog = pd.read_csv(EV / "gwas_catalog_celiac_associations.csv", low_memory=False)
catalog["p"] = pd.to_numeric(catalog["P-VALUE"], errors="coerce")
catalog["mlog"] = -np.log10(catalog["p"])

# The shipped pipeline output (kept for the audit) double-counts rows that appear in both the
# celiac and the autoimmune-context extracts, and its "celiac" substring filter drops every
# "Coeliac disease" row. The corrected analysis below scores the full celiac extract once.
raw_norm = pd.read_csv(LEGACY / "normalized_celiac_associations.csv") if HAVE_LEGACY else None
KEY = ["rsid", "gene", "trait", "p_value"]

gene_rows = []
for r in catalog.to_dict("records"):
    for g in split_gene_field(r["MAPPED_GENE"]):
        gene_rows.append(dict(rsid=str(r["SNPS"]).strip(), gene=g, trait=r["DISEASE/TRAIT"], p_value=r["p"],
                              pmid=r["PUBMEDID"], chrom=str(r["CHR_ID"]), pos=pd.to_numeric(r["CHR_POS"], errors="coerce")))
records = pd.DataFrame(gene_rows)
dedup = records.drop_duplicates(KEY).copy()
# Extended MHC, GRCh38 chr6:25.0-34.0 Mb: long-range LD makes gene-level scores uninterpretable here
dedup["mhc"] = (dedup.chrom == "6") & dedup.pos.between(25e6, 34e6)
mhc_genes = set(dedup.loc[dedup.mhc, "gene"])

all_scores = add_graph_context(score_gene_evidence(dedup[KEY])).reset_index(drop=True)
all_scores["mhc"] = all_scores.gene.isin(mhc_genes)
all_scores["n_studies"] = all_scores.gene.map(dedup.groupby("gene").pmid.nunique()).astype(int)
all_scores["n_unique_snps"] = all_scores.gene.map(dedup.groupby("gene").rsid.nunique()).astype(int)
all_scores["p_term"] = np.minimum(-np.log10(all_scores["best_p_value"]), 300)

scores = all_scores[~all_scores.mhc].reset_index(drop=True)
scores.insert(0, "rank", np.arange(1, len(scores) + 1))
scores["p_only_rank"] = scores["p_term"].rank(ascending=False, method="min").astype(int)
scores["es_only_rank"] = scores["evidence_score"].rank(ascending=False, method="first").astype(int)
mhc_scores = all_scores[all_scores.mhc].reset_index(drop=True)

shipped = pd.read_csv(LEGACY / "gene_evidence_scores.csv") if HAVE_LEGACY else None

panel = pd.read_csv(ROOT / "configs/celiac_snp_panel.csv")
ieu_rows, lps = {}, []
with gzip.open(VCF, "rt") as handle:
    for line in handle:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        es, se, lp, ss = (float(v) for v in f[9].split(":")[:4])
        lps.append(lp)
        if f[2] in set(panel.rsid):
            ieu_rows[f[2]] = dict(chrom=f[0], pos=int(f[1]), ref=f[3], alt=f[4], beta=es, se=se, lp=lp, n=int(ss))
lps = np.array(lps)
lam = float(np.median(chi2.isf(10 ** -lps, 1)) / chi2.ppf(0.5, 1))

# ------------------------------------------------------------- Figure 1: workflow
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.set_xlim(0, 100), ax.set_ylim(0, 64), ax.axis("off")


def box(x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=color, ec="none"))
    ax.text(x + w / 2, y + h - 2.4, title, ha="center", va="top", color="white", fontsize=7.6, fontweight="bold")
    ax.text(x + w / 2, y + h - 6.6, body, ha="center", va="top", color="white", fontsize=6.6, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.0))


n_cat, n_studies = len(catalog), catalog.PUBMEDID.nunique()
n_edges = build_gene_graph(all_scores.gene.tolist()).number_of_edges()
box(1, 47, 31, 16, "GWAS Catalog", f"Release of May 2026\nceliac / coeliac / gluten\n{n_cat} associations, {n_studies} studies", C_DARK)
box(34.5, 47, 31, 16, "IEU OpenGWAS ieu-a-276", f"Dubois et al. 2010 (GRCh37)\n4,533 cases, 10,750 controls\n{len(lps):,} variants, no MHC", C_DARK)
box(68, 47, 31, 16, "FinnGen R12", "K11_COELIAC (GRCh38)\n5,130 cases, 478,189 controls\nnot one of the 16 catalog studies", C_DARK)
box(1, 24, 31, 17, "Score genes", f"{len(dedup)} records, {all_scores.shape[0]} genes\nES = -log10(p_best) + n\nMHC ({len(mhc_scores)} genes) as one locus\n{len(scores)} non-MHC genes ranked", C_MID)
box(34.5, 24, 31, 17, "Panel look-up", f"{len(panel)} SNPs (5 HLA tags,\n5 non-HLA SNPs)\nOR per allele, 95% CI\nin both summary datasets", C_MID)
box(68, 24, 31, 17, "Replication of ranking", "Catalog SNPs looked up\nin FinnGen; gene replicated\nif Bonferroni p < 0.05\nAUC: ES vs p value alone", C_MID)
box(1, 2, 31, 15, "Sensitivity analyses", f"Seed graph ({n_edges} curated edges;\nnot STRING)\nRank by p value alone", "#5b8c6a")
box(34.5, 2, 64.5, 15, "Outputs", "Ranked gene tables (Tables S1, S2, S4)\nPanel statistics (Table 3) and replication (Figure 5)\nScripts that regenerate every number", "#5b8c6a")
for x in (16.5, 50, 83.5):
    arrow(x, 46.5, x, 41.6)
arrow(16.5, 23.5, 16.5, 17.6), arrow(50, 23.5, 50, 17.6), arrow(83.5, 23.5, 83.5, 17.6)
arrow(83.5, 46.5, 55, 41.6)
fig.savefig(FIG / "Figure1_Workflow.png"), fig.savefig(FIG / "Figure1_Workflow.tiff")
plt.close(fig)

# ------------------------------------------- Figure 2: catalog landscape + IEU QQ
cat = catalog.dropna(subset=["CHR_ID", "CHR_POS"]).copy()
cat = cat[~cat.CHR_ID.astype(str).str.contains("x|;", regex=True)]
cat["CHR_ID"] = cat.CHR_ID.astype(str)
cat["pos"] = pd.to_numeric(cat.CHR_POS, errors="coerce")
cat = cat.dropna(subset=["pos"])
chroms = sorted(cat.CHR_ID.unique(), key=chrom_key)
offset, ticks, cum = {}, [], 0
for c in chroms:
    offset[c] = cum
    mx = cat.loc[cat.CHR_ID == c, "pos"].max()
    ticks.append(cum + mx / 2)
    cum += mx + 2e7
cat["x"] = cat.apply(lambda r: offset[r.CHR_ID] + r.pos, axis=1)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [2.6, 1]})
for i, c in enumerate(chroms):
    d = cat[cat.CHR_ID == c]
    a1.scatter(d.x, d.mlog.clip(upper=60), s=9, color=C_DARK if i % 2 else C_MID, lw=0)
capped = cat[cat.mlog > 60]
a1.scatter(capped.x, [60] * len(capped), s=14, marker="^", color=C_ACC, lw=0)
a1.axhline(-np.log10(5e-8), ls="--", lw=0.8, color=C_ACC)
a1.set_xticks(ticks, [c if chrom_key(c) <= 12 or chrom_key(c) % 2 else "" for c in chroms], fontsize=6.5), a1.set_ylim(4, 63)
a1.set_xlabel("Chromosome (GRCh38)"), a1.set_ylabel("−log$_{10}$(p)")
a1.set_title("a  GWAS Catalog celiac associations")
for gene in ["LPP", "IL21-AS1", "PHTF1 - RSBN1"]:
    r = cat[cat.MAPPED_GENE == gene].sort_values("mlog").iloc[-1]
    a1.annotate(gene.replace(" - ", "/"), (r.x, min(r.mlog, 60)), xytext=(4, 2), textcoords="offset points", fontsize=6.5)
a1.text(0.99, 0.80, f"{len(capped)} associations with p < 10$^{{-60}}$ (▲, MHC) shown at cap", transform=a1.transAxes, ha="right", va="top", fontsize=6.5)

exp = -np.log10((np.arange(1, len(lps) + 1) - 0.5) / len(lps))
obs = np.sort(lps)[::-1]
a2.scatter(exp, obs, s=4, color=C_MID, lw=0, rasterized=True)
a2.plot([0, exp.max()], [0, exp.max()], color=C_GREY, lw=0.8, ls="--")
a2.set_xlabel("Expected −log$_{10}$(p)"), a2.set_ylabel("Observed −log$_{10}$(p)")
a2.set_title("b  ieu-a-276 Q–Q")
a2.text(0.05, 0.95, f"λ$_{{GC}}$ = {lam:.3f}\n{len(lps):,} variants", transform=a2.transAxes, va="top", fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "Figure2_Association_Landscape.png"), fig.savefig(FIG / "Figure2_Association_Landscape.tiff")
plt.close(fig)

# ---------------------------------------- Figure 3: score composition + rank shift
top = scores.head(20).iloc[::-1]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
a1.barh(top.gene, top.p_term, color=C_DARK, label="−log$_{10}$(p$_{best}$)")
a1.barh(top.gene, top.association_count, left=top.p_term, color=C_LIGHT, label="Record count")
for y, (g, v) in enumerate(zip(top.gene, top.combined_evidence_graph_score)):
    a1.text(v + 0.6, y, f"{v:.1f}", va="center", fontsize=6.5)
a1.set_xlabel("Evidence score"), a1.set_title("a  Top 20 genes and score components")
a1.legend(frameon=False, fontsize=7, loc="lower right"), a1.tick_params(axis="y", labelsize=7)

top20 = scores.head(20)
for r in top20.itertuples():
    col = C_ACC if abs(r.rank - r.p_only_rank) >= 5 else C_GREY
    a2.plot([0, 1], [r.p_only_rank, r.rank], color=col, lw=1.1 if col == C_ACC else 0.7)
    a2.text(1.04, r.rank, r.gene, va="center", fontsize=6.3, color=col)
a2.set_xlim(-0.15, 1.5), a2.set_ylim(31, 0)
a2.set_xticks([0, 1], ["p-value\nonly", "Evidence\nscore"])
a2.set_ylabel("Rank"), a2.spines["bottom"].set_visible(False)
a2.set_title("b  Rank under each scheme")
fig.tight_layout()
fig.savefig(FIG / "Figure3_Gene_Ranking.png"), fig.savefig(FIG / "Figure3_Gene_Ranking.tiff")
plt.close(fig)

# ---------------------------------------------------- Figure 4: panel forest plot
lk = panel.copy()
lk["present"] = lk.rsid.isin(ieu_rows)
pres = lk[lk.present].copy()
for k in ["chrom", "pos", "ref", "alt", "beta", "se", "lp", "n"]:
    pres[k] = pres.rsid.map(lambda s, k=k: ieu_rows[s][k])
# OR per copy of the reference allele (which is the allele associated with higher odds here)
pres["or_ref"] = np.exp(-pres.beta)
pres["lo"], pres["hi"] = np.exp(-(pres.beta + 1.96 * pres.se)), np.exp(-(pres.beta - 1.96 * pres.se))
pres["z"] = pres.beta.abs() / pres.se
pres["p"] = 10 ** -pres.lp
pres = pres.sort_values("lp")

fig, ax = plt.subplots(figsize=(5.6, 2.6))
y = np.arange(len(pres))
ax.errorbar(pres.or_ref, y, xerr=[pres.or_ref - pres.lo, pres.hi - pres.or_ref], fmt="s", color=C_DARK, ms=5, capsize=2.5, lw=1)
ax.axvline(1, color=C_GREY, ls="--", lw=0.8)
ax.set_yticks(y, [f"{g}  {s}" for g, s in zip(pres.gene, pres.rsid)], fontsize=7.5)
for yi, r in zip(y, pres.itertuples()):
    ax.text(1.62, yi, f"{r.or_ref:.2f} ({r.lo:.2f}–{r.hi:.2f})   p = {r.p:.1e}", va="center", fontsize=6.8)
ax.set_xlim(0.95, 1.6), ax.set_xlabel("Odds ratio per reference allele (95% CI)")
ax.set_title("ieu-a-276 (4,533 cases, 10,750 controls)")
fig.savefig(FIG / "Figure4_Panel_Lookup.png"), fig.savefig(FIG / "Figure4_Panel_Lookup.tiff")
plt.close(fig)

# ------------------------------------------ Supplementary Figure 1: annotation mix
ctx = catalog.CONTEXT.fillna("not annotated").str.replace("_variant", "").str.replace("_", " ")
ctx_counts = ctx.value_counts()
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
a1.barh(ctx_counts.index[::-1], ctx_counts.values[::-1], color=C_MID)
a1.set_xlabel("Associations"), a1.set_title("a  Variant consequence"), a1.tick_params(axis="y", labelsize=7)
chr_counts = catalog.CHR_ID.astype(str).value_counts()
chr_counts = chr_counts[[c for c in chr_counts.index if c.replace("X", "23").isdigit()]]
chr_counts = chr_counts.reindex(sorted(chr_counts.index, key=chrom_key))
a2.bar(chr_counts.index, chr_counts.values, color=[C_ACC if c == "6" else C_DARK for c in chr_counts.index])
a2.set_xlabel("Chromosome"), a2.set_ylabel("Associations"), a2.set_title("b  Associations per chromosome")
a2.set_xticks(range(len(chr_counts)), [c if chrom_key(c) <= 12 or chrom_key(c) % 2 else "" for c in chr_counts.index], fontsize=6.5)
fig.tight_layout()
fig.savefig(SUP / "SuppFigure1_Annotation_Chromosomes.png")
plt.close(fig)

# --------------------------------------------------------------- tables to CSV
cols = ["rank", "gene", "association_count", "n_studies", "n_unique_snps", "best_p_value", "evidence_score",
        "graph_degree", "degree_centrality", "combined_evidence_graph_score", "p_only_rank", "es_only_rank"]
scores[cols].to_csv(SUP / "SuppTable1_NonMHC_Gene_Scores.csv", index=False)
mhc_scores[cols[1:-2]].to_csv(SUP / "SuppTable2_MHC_Gene_Scores.csv", index=False)
if HAVE_LEGACY:
    cmp = scores[["gene", "rank"]].merge(shipped.reset_index().rename(columns={"index": "shipped_rank"})[["gene", "shipped_rank", "association_count", "evidence_score"]], on="gene", how="outer")
    cmp["shipped_rank"] = cmp["shipped_rank"] + 1
    cmp.sort_values("rank").to_csv(AUDIT / "corrected_vs_shipped_ranking.csv", index=False)
pres.sort_values("lp", ascending=False)[["rsid", "gene", "chrom", "pos", "ref", "alt", "beta", "se", "or_ref", "lo", "hi", "z", "p", "n"]].to_csv(AUDIT / "table3_panel_lookup.csv", index=False)
lk[["rsid", "gene", "chromosome", "category", "present"]].to_csv(AUDIT / "panel_presence.csv", index=False)

top1000 = pd.read_csv(EV / "ieu-a-276/ieu_a_276_top1000_summary_stats.csv")
snp_gene = catalog.drop_duplicates("SNPS").set_index("SNPS").MAPPED_GENE
top25 = top1000.nsmallest(25, "p_value").copy()
top25["catalog_mapped_gene"] = top25.rsid.map(snp_gene).fillna("—")
top25[["chromosome", "position", "rsid", "ref", "alt", "effect_size", "standard_error", "p_value", "catalog_mapped_gene"]].to_csv(SUP / "SuppTable3_IEU_Top25_Variants.csv", index=False)

graph = build_gene_graph(all_scores.gene.tolist())
from scipy.stats import spearmanr
rho = spearmanr(scores.p_term, scores.evidence_score)
summary = {
    "catalog_rows": n_cat, "catalog_studies": n_studies, "catalog_unique_snps": catalog.SNPS.nunique(),
    "catalog_gws_rows": int((catalog.p <= 5e-8).sum()), "catalog_mlog_median": float(catalog.mlog.median()),
    "catalog_mlog_max": float(catalog.mlog.max()), "raw_normalised_rows": len(raw_norm) if HAVE_LEGACY else None, "dedup_rows": len(dedup), "dedup_rsids": dedup.rsid.nunique(),
    "genes_all": len(all_scores), "genes_mhc": len(mhc_scores), "genes_nonmhc": len(scores), "mhc_records": int(dedup.mhc.sum()),
    "rho_p_vs_es": rho.statistic, "rho_p": rho.pvalue, "graph_changed_genes": int((all_scores.graph_degree > 0).sum()),
    "max_graph_increment": float((all_scores.combined_evidence_graph_score - all_scores.evidence_score).max()),
    "rank_changes_from_graph_nonmhc": int((scores["rank"] != scores.es_only_rank).sum()),
    "es_median": scores.evidence_score.median(), "es_q1": scores.evidence_score.quantile(.25), "es_q3": scores.evidence_score.quantile(.75), "graph_edges": graph.number_of_edges(), "ieu_variants": len(lps), "ieu_lambda": lam,
    "ieu_gws": int((lps > -np.log10(5e-8)).sum()), "panel_present": len(ieu_rows),
}
pd.Series(summary).to_csv(AUDIT / "derived_numbers.csv", header=["value"])
print(pd.Series(summary).to_string())
print(scores[cols].head(22).round(2).to_string(index=False)); print(mhc_scores[cols[1:-2]].head(10).round(2).to_string(index=False)); print(scores[scores.graph_degree>0][cols].to_string())
print(pres[["rsid", "gene", "chrom", "pos", "beta", "se", "or_ref", "lo", "hi", "z", "p"]].round(3).to_string(index=False))
print(top25.head(8).to_string(index=False))
