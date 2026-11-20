"""Step 4: leave-one-donor-out (LODO) robustness of pseudobulk DEGs.

For each cell type with enough DEGs, PyDESeq2 is re-run once per donor with that donor removed.
A DEG is "retained" in a run if it stays significant (same FDR and |log2FC| thresholds) with the
same direction. DEGs retained in every run do not depend on any single donor.
"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _de_module():
    spec = importlib.util.spec_from_file_location("pseudobulk_de", Path(__file__).with_name("02_pseudobulk_de.py"))
    if not spec.origin or not Path(spec.origin).exists():  # Snakemake copies scripts; fall back to the repo path
        spec = importlib.util.spec_from_file_location("pseudobulk_de", "workflow/scripts/02_pseudobulk_de.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def is_deg(df: pd.DataFrame, padj: float, lfc: float) -> pd.Series:
    return (df["padj"] < padj) & (df["log2FoldChange"].abs() >= lfc)


def leave_one_donor_out(counts: pd.DataFrame, meta: pd.DataFrame, full_de: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    m2 = _de_module()
    pb, rob = cfg["pseudobulk"], cfg["robustness"]
    padj, lfc = pb["padj_threshold"], pb["lfc_threshold"]
    gene_rows, summary = [], []
    for ct, full in full_de.groupby("cell_type"):
        degs = full[is_deg(full, padj, lfc)].set_index("gene")
        if len(degs) < int(rob["min_degs"]):
            continue
        m_ct = meta[meta["cell_type"] == ct]
        retained = pd.Series(0, index=degs.index)
        runs = 0
        for donor in sorted(m_ct["donor"].unique()):
            m = m_ct[m_ct["donor"] != donor]
            if m["condition"].value_counts().min() < int(pb["min_donors_per_group"]) - 1:
                continue
            c = m2.filter_genes(counts.loc[m.index], m, int(pb["min_gene_count"]))
            res = m2.run_deseq2(c, m, int(cfg.get("threads", 1))).set_index("gene")
            res = res.reindex(degs.index)
            ok = is_deg(res.fillna({"padj": 1.0, "log2FoldChange": 0.0}), padj, lfc) & (
                np.sign(res["log2FoldChange"].fillna(0)) == np.sign(degs["log2FoldChange"]))
            retained += ok.astype(int)
            runs += 1
        frac = retained / runs if runs else retained * np.nan
        gene_rows.append(pd.DataFrame({"cell_type": ct, "gene": degs.index, "log2FoldChange": degs["log2FoldChange"].to_numpy(),
                                       "padj": degs["padj"].to_numpy(), "n_runs": runs, "n_retained": retained.to_numpy(),
                                       "fraction_retained": frac.to_numpy()}))
        summary.append({"cell_type": ct, "n_degs": len(degs), "n_runs": runs,
                        "retained_all_runs": int((frac == 1).sum()), "pct_retained_all_runs": 100 * float((frac == 1).mean()),
                        "retained_90pct_runs": int((frac >= 0.9).sum()), "median_fraction_retained": float(frac.median())})
    genes = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    return genes, pd.DataFrame(summary)


def plot_fig4(genes: pd.DataFrame, summary: pd.DataFrame, out_path: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), dpi=180)
    if genes.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No cell type met the DEG threshold for LODO", ha="center", va="center")
    else:
        ax = axes[0]
        for ct, sub in genes.groupby("cell_type"):
            x = np.sort(sub["fraction_retained"].to_numpy())[::-1]
            ax.plot(np.arange(1, len(x) + 1) / len(x) * 100, x * 100, lw=2, label=ct.replace(" of endometrium", ""))
        ax.set_xlabel("DEGs ranked by retention (%)")
        ax.set_ylabel("Leave-one-donor-out runs in which the DEG is retained (%)")
        ax.set_title("Robustness of DEGs to removing any single donor", fontweight="bold", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax = axes[1]
        s = summary.copy()
        y = np.arange(len(s))
        ax.barh(y, s["n_degs"], color="#cbd5e0", label="DEGs (all donors)")
        ax.barh(y, s["retained_90pct_runs"], color="#718096", label="retained in ≥ 90% of runs")
        ax.barh(y, s["retained_all_runs"], color="#2d3748", label="retained in every run")
        ax.set_yticks(y)
        ax.set_yticklabels([c.replace(" of endometrium", "") for c in s["cell_type"]], fontsize=8)
        ax.set_xlabel("Number of DEGs")
        ax.set_title("DEGs that do not depend on any single donor", fontweight="bold", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_step(pb_counts, pb_meta, de_path, cfg, out_genes, out_summary, out_fig) -> None:
    counts = pd.read_csv(pb_counts, index_col=0)
    meta = pd.read_csv(pb_meta, index_col=0)
    full = pd.read_csv(de_path)
    genes, summary = leave_one_donor_out(counts, meta, full, cfg)
    for df, p in ((genes, out_genes), (summary, out_summary)):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)
    plot_fig4(genes, summary, out_fig)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_step(snakemake.input.pb_counts, snakemake.input.pb_meta, snakemake.input.de, snakemake.params.cfg,
                 snakemake.output.genes, snakemake.output.summary, snakemake.output.fig4)
