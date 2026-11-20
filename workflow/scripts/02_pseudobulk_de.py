"""Step 2: donor-level pseudobulk differential expression per cell type with PyDESeq2.

Each sample is one donor x cell type (raw counts summed over all nuclei). The design is
`~ condition` (programmed vs natural cycle), with donors as the unit of replication; this avoids
the inflated false-positive rate of treating individual nuclei as independent replicates.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def filter_genes(counts: pd.DataFrame, meta: pd.DataFrame, min_count: int) -> pd.DataFrame:
    """Keep genes with >= min_count reads in at least as many samples as the smaller group."""
    n_min = int(meta["condition"].value_counts().min())
    keep = (counts >= min_count).sum(axis=0) >= n_min
    return counts.loc[:, keep]


def run_deseq2(counts: pd.DataFrame, meta: pd.DataFrame, n_cpus: int = 1) -> pd.DataFrame:
    """PyDESeq2 Wald test for condition 'test' vs 'reference'."""
    inference = DefaultInference(n_cpus=n_cpus)
    md = meta[["condition"]].copy()
    try:
        dds = DeseqDataSet(counts=counts, metadata=md, design="~ condition", refit_cooks=True,
                           inference=inference, quiet=True)
    except TypeError:  # PyDESeq2 < 0.5
        dds = DeseqDataSet(counts=counts, metadata=md, design_factors="condition", refit_cooks=True,
                           inference=inference, quiet=True)
    dds.deseq2()
    res = DeseqStats(dds, contrast=["condition", "test", "reference"], inference=inference, quiet=True)
    res.summary()
    return res.results_df.reset_index().rename(columns={"index": "gene"})


def run_all(counts: pd.DataFrame, meta: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    pb = cfg["pseudobulk"]
    de_frames, summary = [], []
    for ct, m in meta.groupby("cell_type"):
        n_test = int((m["condition"] == "test").sum())
        n_ref = int((m["condition"] == "reference").sum())
        row = {"cell_type": ct, "n_test_donors": n_test, "n_reference_donors": n_ref,
               "n_genes_tested": 0, "n_up": 0, "n_down": 0, "tested": False}
        if min(n_test, n_ref) < int(pb["min_donors_per_group"]):
            summary.append(row)
            continue
        c = filter_genes(counts.loc[m.index], m, int(pb["min_gene_count"]))
        df = run_deseq2(c, m, int(cfg.get("threads", 1)))
        df.insert(0, "cell_type", ct)
        df["n_test_donors"], df["n_reference_donors"] = n_test, n_ref
        sig = df[(df["padj"] < pb["padj_threshold"]) & (df["log2FoldChange"].abs() >= pb["lfc_threshold"])]
        row.update({"n_genes_tested": len(df), "n_up": int((sig["log2FoldChange"] > 0).sum()),
                    "n_down": int((sig["log2FoldChange"] < 0).sum()), "tested": True})
        summary.append(row)
        de_frames.append(df)
    de = pd.concat(de_frames, ignore_index=True) if de_frames else pd.DataFrame()
    return de, pd.DataFrame(summary)


def plot_fig2(de: pd.DataFrame, summary: pd.DataFrame, cfg: dict, out_path: str) -> None:
    pb = cfg["pseudobulk"]
    d = cfg["dataset"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), dpi=180, gridspec_kw={"width_ratios": [1, 1.2]})
    ax = axes[0]
    s = summary[summary["tested"]].sort_values("n_up", ascending=True)
    y = np.arange(len(s))
    ax.barh(y, s["n_up"], color="#c05621", label=f"Higher in {d['test_level']}")
    ax.barh(y, -s["n_down"], color="#2b6cb0", label=f"Higher in {d['reference_level']}")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({a} vs {b} donors)" for c, a, b in zip(s["cell_type"], s["n_test_donors"], s["n_reference_donors"])], fontsize=7.5)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"DEGs (FDR < {pb['padj_threshold']}, |log2FC| ≥ {pb['lfc_threshold']})")
    ax.set_title("Pseudobulk DEGs per cell type", fontweight="bold", fontsize=10)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    if not de.empty and s["n_up"].add(s["n_down"]).max() > 0:
        top_ct = s.assign(n=s["n_up"] + s["n_down"]).sort_values("n").iloc[-1]["cell_type"]
        v = de[de["cell_type"] == top_ct].dropna(subset=["padj"]).copy()
        v["nlp"] = -np.log10(v["padj"].clip(lower=1e-300))
        sig = (v["padj"] < pb["padj_threshold"]) & (v["log2FoldChange"].abs() >= pb["lfc_threshold"])
        col = np.where(~sig, "#cbd5e0", np.where(v["log2FoldChange"] > 0, "#c05621", "#2b6cb0"))
        ax.scatter(v["log2FoldChange"], v["nlp"], s=5, c=col, rasterized=True)
        for _, r in v[sig].nlargest(12, "nlp").iterrows():
            ax.annotate(r["gene"], (r["log2FoldChange"], r["nlp"]), fontsize=7, xytext=(3, 2), textcoords="offset points")
        ax.axhline(-np.log10(pb["padj_threshold"]), color="#a0aec0", ls="--", lw=0.8)
        ax.set_xlabel(f"log2 fold change ({d['test_level']} / {d['reference_level']})")
        ax.set_ylabel("-log10 adjusted p-value")
        ax.set_title(f"{top_ct}", fontweight="bold", fontsize=10)
    else:
        ax.text(0.5, 0.5, "No DEGs at the chosen thresholds", ha="center", va="center")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_step(pb_counts: str, pb_meta: str, cfg: dict, out_de: str, out_summary: str, out_fig: str) -> None:
    counts = pd.read_csv(pb_counts, index_col=0)
    meta = pd.read_csv(pb_meta, index_col=0)
    de, summary = run_all(counts, meta, cfg)
    for df, p in ((de, out_de), (summary, out_summary)):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)
    plot_fig2(de, summary, cfg, out_fig)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_step(snakemake.input.pb_counts, snakemake.input.pb_meta, snakemake.params.cfg,
                 snakemake.output.de, snakemake.output.summary, snakemake.output.fig2)
