"""Independent replication in FinnGen R12 (K11_COELIAC: 5,130 cases, 478,189 controls).

Queries the public FinnGen summary statistics remotely with tabix (no full download) for
(a) the 10 candidate-panel SNPs and (b) every GWAS Catalog celiac SNP with a GRCh38 position,
then asks whether the evidence score predicts FinnGen replication better than p value alone.
Run from the project root after make_figures.py.
"""
from __future__ import annotations

import io
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deepg2p_celiac.analysis.evidence import split_gene_field  # noqa: E402

OUT = ROOT / "results"
(OUT / "audit").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
FINNGEN = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_K11_COELIAC.gz"
CACHE = OUT / "audit/finngen_hits.tsv"
rng = np.random.default_rng(2026)

catalog = pd.read_csv(ROOT / "data/processed/open_source_evidence/gwas_catalog_celiac_associations.csv", low_memory=False)
catalog["pos"] = pd.to_numeric(catalog.CHR_POS, errors="coerce")
panel = pd.read_csv(ROOT / "configs/celiac_snp_panel.csv")

# GRCh38 positions for panel SNPs (Ensembl REST); catalog positions are already GRCh38
req = urllib.request.Request("https://rest.ensembl.org/variation/homo_sapiens", method="POST",
                             data=json.dumps({"ids": panel.rsid.tolist()}).encode(),
                             headers={"Content-Type": "application/json", "Accept": "application/json"})
ens = json.load(urllib.request.urlopen(req, timeout=60))
panel_pos = {k: next((m["seq_region_name"], int(m["start"])) for m in v["mappings"]
                     if m["assembly_name"] == "GRCh38" and len(m["seq_region_name"]) < 3) for k, v in ens.items()}

cat_pos = catalog.dropna(subset=["pos"])[["SNPS", "CHR_ID", "pos"]].drop_duplicates("SNPS")
targets = {r.SNPS: (str(r.CHR_ID), int(r.pos)) for r in cat_pos.itertuples()} | panel_pos

if not CACHE.exists():
    regions = "\n".join(f"{c.replace('X', '23')}\t{p}\t{p}" for c, p in sorted(set(targets.values()))) + "\n"  # FinnGen contigs: 1..22, 23
    reg_file = OUT / "audit/finngen_regions.tsv"
    reg_file.write_text(regions)
    res = subprocess.run(["tabix", "-R", str(reg_file), FINNGEN], capture_output=True, text=True, check=True)
    header = "chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls\n"
    CACHE.write_text(header + res.stdout)
fg = pd.read_csv(CACHE, sep="\t")
fg = fg.drop_duplicates()


def lookup(rsid):
    """Match by rsID first (handles multi-allelic sites), else by position if only one record."""
    hit = fg[fg.rsids.fillna("").str.split(",").apply(lambda ids: rsid in ids)]
    if hit.empty and rsid in targets:
        c, p = targets[rsid]
        hit = fg[(fg.chrom.astype(str) == c.replace("X", "23")) & (fg.pos == p)]
        hit = hit if len(hit) == 1 else hit.iloc[0:0]
    return hit.iloc[0] if len(hit) else None


# ---------------------------------------------------------------- (a) panel SNPs
ieu = pd.read_csv(OUT / "audit/table3_panel_lookup.csv").set_index("rsid")
rows = []
for r in panel.itertuples():
    h = lookup(r.rsid)
    if h is None:
        rows.append(dict(rsid=r.rsid, gene=r.gene, present=False)); continue
    beta, se = float(h.beta), float(h.sebeta)
    rec = dict(rsid=r.rsid, gene=r.gene, present=True, chrom=h.chrom, pos38=int(h.pos), ref=h.ref, alt=h.alt,
               beta=beta, se=se, p=float(h.pval), mlogp=float(h.mlogp), af_alt=float(h.af_alt))
    if r.rsid in ieu.index:  # harmonize to the ieu-a-276 reference allele (higher-odds allele there)
        i = ieu.loc[r.rsid]
        if h.ref == i.ref and h.alt == i.alt:
            sign = -1
        elif h.ref == i.alt and h.alt == i.ref:
            sign = 1
        else:
            sign = np.nan
        rec.update(ieu_ref=i.ref, harmonized_log_or=sign * beta,
                   or_ieu_ref=np.exp(sign * beta), lo=np.exp(sign * beta - 1.96 * se), hi=np.exp(sign * beta + 1.96 * se),
                   same_direction=bool(sign * beta > 0) if sign == sign else None)
    else:  # HLA tags: report per FinnGen alt allele
        rec.update(or_alt=np.exp(beta), lo=np.exp(beta - 1.96 * se), hi=np.exp(beta + 1.96 * se))
    rows.append(rec)
panel_rep = pd.DataFrame(rows)
panel_rep.to_csv(OUT / "audit/finngen_panel_replication.csv", index=False)

# ------------------------------------------------- (b) catalog SNPs and gene ranking
snp_rep = []
for rsid in cat_pos.SNPS:
    h = lookup(rsid)
    snp_rep.append(dict(rsid=rsid, fg_p=float(h.pval) if h is not None else np.nan,
                        fg_mlogp=float(h.mlogp) if h is not None else np.nan))
snp_rep = pd.DataFrame(snp_rep)

scores = pd.read_csv(OUT / "supplementary/SuppTable1_NonMHC_Gene_Scores.csv")
gene_snps = {}
for r in catalog.to_dict("records"):
    for g in split_gene_field(r["MAPPED_GENE"]):
        gene_snps.setdefault(g, set()).add(str(r["SNPS"]))
mlog = snp_rep.set_index("rsid").fg_mlogp


def gene_rep(g):
    vals = [mlog.get(s, np.nan) for s in gene_snps.get(g, [])]
    vals = [v for v in vals if v == v]
    if not vals:
        return pd.Series(dict(fg_snps=0, fg_best_mlogp=np.nan, fg_bonf_p=np.nan))
    best = max(vals)
    return pd.Series(dict(fg_snps=len(vals), fg_best_mlogp=best, fg_bonf_p=min(1.0, 10 ** -best * len(vals))))


scores = pd.concat([scores, scores.gene.apply(gene_rep)], axis=1)
scores["replicated"] = scores.fg_bonf_p < 0.05
scores.to_csv(OUT / "audit/finngen_gene_replication.csv", index=False)
ev = scores.dropna(subset=["fg_best_mlogp"]).copy()
ev["p_term"] = -np.log10(ev.best_p_value)


def auc(score, label):
    pos, neg = score[label], score[~label]
    return mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


# Genes that share a catalog SNP (e.g. the two flanking genes of an intergenic SNP) are not
# independent observations, so the bootstrap resamples connected clusters of genes.
import networkx as nx
G = nx.Graph()
G.add_nodes_from(ev.gene)
for g in ev.gene:
    for s_ in gene_snps.get(g, []):
        G.add_edge(g, "snp:" + s_)
cluster_of = {g: i for i, comp in enumerate(nx.connected_components(G)) for g in comp if not g.startswith("snp:")}
ev["cluster"] = ev.gene.map(cluster_of)
clusters = ev.groupby("cluster").groups


def boot_diff(n=2000):
    diffs = []
    keys = np.array(list(clusters))
    for _ in range(n):
        pick = rng.choice(keys, len(keys))
        b = ev.loc[np.concatenate([clusters[k] for k in pick])]
        if b.replicated.nunique() < 2:
            continue
        diffs.append(auc(b.evidence_score, b.replicated) - auc(b.p_term, b.replicated))
    return np.percentile(diffs, [2.5, 97.5])


top20 = ev.nsmallest(20, "rank")
res = {
    "finngen_cases": 5130, "finngen_controls": 478189,
    "catalog_snps_with_pos": len(cat_pos), "catalog_snps_found": int(snp_rep.fg_p.notna().sum()),
    "catalog_snps_p05": int((snp_rep.fg_p < 0.05).sum()), "catalog_snps_gws": int((snp_rep.fg_p < 5e-8).sum()),
    "nonmhc_genes_evaluable": len(ev), "nonmhc_genes_replicated": int(ev.replicated.sum()),
    "top20_replicated": int(top20.replicated.sum()), "top20_evaluable": len(top20),
    "rest_replicated": int(ev[ev["rank"] > 20].replicated.sum()), "rest_evaluable": int((ev["rank"] > 20).sum()),
    "rho_es_vs_fg": spearmanr(ev.evidence_score, ev.fg_best_mlogp).statistic,
    "rho_es_vs_fg_p": spearmanr(ev.evidence_score, ev.fg_best_mlogp).pvalue,
    "rho_pterm_vs_fg": spearmanr(ev.p_term, ev.fg_best_mlogp).statistic,
    "rho_pterm_vs_fg_p": spearmanr(ev.p_term, ev.fg_best_mlogp).pvalue,
    "rho_count_vs_fg": spearmanr(ev.association_count, ev.fg_best_mlogp).statistic,
    "rho_count_vs_fg_p": spearmanr(ev.association_count, ev.fg_best_mlogp).pvalue,
    "auc_es": auc(ev.evidence_score, ev.replicated), "auc_pterm": auc(ev.p_term, ev.replicated),
    "auc_count": auc(ev.association_count, ev.replicated),
}
res["auc_diff_ci_lo"], res["auc_diff_ci_hi"] = boot_diff()
res["n_clusters"] = len(clusters)
pd.Series(res).to_csv(OUT / "audit/finngen_summary.csv", header=["value"])
pd.set_option("display.width", 220)
print(pd.Series(res).to_string())
print(panel_rep.round(4).to_string())
print(ev[["rank", "gene", "evidence_score", "fg_snps", "fg_best_mlogp", "fg_bonf_p", "replicated"]].head(20).round(3).to_string())

# ------------------------------------------- (c) direction check for catalog SNPs
COMP = str.maketrans("ACGT", "TGCA")
dir_rows = []
for r in catalog.to_dict("records"):
    ra = str(r["STRONGEST SNP-RISK ALLELE"]).split("-")[-1].strip()
    or_beta = pd.to_numeric(r["OR or BETA"], errors="coerce")
    ci = str(r["95% CI (TEXT)"])
    h = lookup(str(r["SNPS"]))
    # only odds ratios (no "increase/decrease" beta text), a usable risk allele, a non-palindromic FinnGen site
    if h is None or ra not in list("ACGT") or not (or_beta == or_beta) or "increase" in ci or "decrease" in ci or or_beta <= 0:
        continue
    if {h.ref, h.alt} in ({"A", "T"}, {"C", "G"}) or len(h.ref) != 1 or len(h.alt) != 1:
        continue
    allele = ra if ra in (h.ref, h.alt) else ra.translate(COMP) if ra.translate(COMP) in (h.ref, h.alt) else None
    if allele is None:
        continue
    catalog_up = or_beta > 1                     # risk allele raises odds in the discovery study
    fg_up = (h.beta > 0) == (allele == h.alt)    # same allele raises odds in FinnGen
    dir_rows.append(dict(rsid=r["SNPS"], study=r["FIRST AUTHOR"], fg_p=float(h.pval), concordant=catalog_up == fg_up))
dr = pd.DataFrame(dir_rows).drop_duplicates("rsid")
sig = dr[dr.fg_p < 0.05]
from scipy.stats import binomtest
bt = binomtest(int(sig.concordant.sum()), len(sig), 0.5)
dsum = {"dir_snps_checked": len(dr), "dir_snps_p05": len(sig), "dir_concordant_p05": int(sig.concordant.sum()),
        "dir_binom_p": bt.pvalue, "dir_all_concordant": int(dr.concordant.sum())}
print(pd.Series(dsum).to_string())
print(sig[~sig.concordant])
s = pd.read_csv(OUT / "audit/finngen_summary.csv", index_col=0)["value"]
pd.concat([s, pd.Series(dsum)]).to_csv(OUT / "audit/finngen_summary.csv", header=["value"])
dr.to_csv(OUT / "audit/finngen_direction_check.csv", index=False)


# ----------------------------------------------------------- Figure 5: FinnGen replication
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 10, "axes.titlelocation": "left",
                     "savefig.dpi": 300, "savefig.bbox": "tight"})
C_DARK, C_MID, C_ACC, C_GREY = "#1f3b5a", "#2a7f9e", "#c8553d", "#8c8c8c"
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1, 1.15]})
pr = panel_rep[panel_rep.or_ieu_ref.notna()].merge(ieu[["or_ref", "lo", "hi"]].rename(columns={"lo": "ilo", "hi": "ihi"}),
                                                    left_on="rsid", right_index=True).sort_values("or_ref")
y = np.arange(len(pr))
a1.errorbar(pr.or_ref, y + 0.15, xerr=[pr.or_ref - pr.ilo, pr.ihi - pr.or_ref], fmt="s", ms=4, color=C_DARK, capsize=2, lw=1, label="ieu-a-276")
a1.errorbar(pr.or_ieu_ref, y - 0.15, xerr=[pr.or_ieu_ref - pr.lo, pr.hi - pr.or_ieu_ref], fmt="o", ms=4, color=C_ACC, capsize=2, lw=1, label="FinnGen R12")
a1.axvline(1, color=C_GREY, ls="--", lw=0.8)
a1.set_yticks(y, [f"{g}\n{r}" for g, r in zip(pr.gene, pr.rsid)], fontsize=7)
a1.set_xlabel("OR per ieu-a-276 reference allele (95% CI)")
a1.set_title("a  Non-HLA panel SNPs")
a1.legend(frameon=False, fontsize=7, loc="lower right")
col = np.where(ev.replicated, C_MID, C_ACC)
a2.scatter(ev.evidence_score, ev.fg_best_mlogp, s=12, c=col, lw=0)
a2.axhline(-np.log10(0.05), color=C_GREY, ls=":", lw=0.8)
a2.set_xlabel("Evidence score (non-MHC genes)")
a2.set_ylabel("Best FinnGen −log$_{10}$(p)")
a2.set_title("b  Gene ranking vs FinnGen")
for g in ["LPP", "RSBN1", "CCR3", "TAGAP", "IL21-AS1"]:
    r = ev[ev.gene == g].iloc[0]
    a2.annotate(g, (r.evidence_score, r.fg_best_mlogp), xytext=(3, 2), textcoords="offset points", fontsize=6.5)
a2.text(0.03, 0.97, f"Spearman ρ = {res['rho_es_vs_fg']:.2f}\nreplicated (blue): {res['nonmhc_genes_replicated']}/{res['nonmhc_genes_evaluable']}",
        transform=a2.transAxes, ha="left", va="top", fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "figures/Figure5_FinnGen_Replication.png"), fig.savefig(OUT / "figures/Figure5_FinnGen_Replication.tiff")
