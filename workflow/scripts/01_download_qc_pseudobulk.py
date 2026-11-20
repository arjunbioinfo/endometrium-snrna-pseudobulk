"""Step 1: download the public endometrium snRNA-seq dataset, run QC, summarise composition,
integrate a donor-balanced subsample with Harmony for visualisation, and aggregate raw counts
into donor x cell-type pseudobulk samples (using all nuclei).

Dataset: Huang et al., JCI Insight 2026 (doi:10.1172/jci.insight.202443), window-of-implantation
endometrium from IVF patients prepared with a programmed (hormone-replacement) or natural cycle;
single-nucleus RNA-seq (10x Flex), 16 donors, distributed under CC BY 4.0 by CZ CELLxGENE Discover.
"""

import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# --------------------------------------------------------------------------------------
# Data access
# --------------------------------------------------------------------------------------
def download_dataset(url: str, dest: Path, timeout: int = 600) -> Path:
    """Download the CELLxGENE h5ad once (atomic rename) and reuse the cached copy."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "endometrium-snrna-pseudobulk/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 22)
    with open(tmp, "rb") as fh:
        if fh.read(8) != b"\x89HDF\r\n\x1a\n":
            tmp.unlink()
            raise RuntimeError(f"{url} did not return an HDF5/.h5ad file")
    tmp.replace(dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_counts(path: Path, cfg: dict) -> ad.AnnData:
    """Return an AnnData with raw integer counts in .X, gene symbols as var_names and
    standardised obs columns: donor, condition, cell_type."""
    d = cfg["dataset"]
    full = ad.read_h5ad(path)
    counts = full.raw.X if full.raw is not None else full.X
    var = (full.raw.var if full.raw is not None else full.var).copy()
    symbols = var[d.get("symbol_column", "feature_name")].astype(str).to_numpy()
    obs = full.obs[[d["donor_column"], d["condition_column"], d["celltype_column"]]].copy()
    obs.columns = ["donor", "condition", "cell_type"]
    for c in obs.columns:
        obs[c] = obs[c].astype(str)
    adata = ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32), obs=obs,
                       var=pd.DataFrame(index=pd.Index(symbols, name="gene")))
    adata.var_names_make_unique()
    adata.obs["condition"] = adata.obs["condition"].map(
        {d["test_level"]: "test", d["reference_level"]: "reference"}
    )
    if adata.obs["condition"].isna().any():
        raise ValueError("Unexpected condition levels; check dataset.test_level / reference_level")
    adata.obs["condition_label"] = np.where(adata.obs["condition"] == "test", d["test_level"], d["reference_level"])
    return adata


# --------------------------------------------------------------------------------------
# QC, composition and pseudobulk
# --------------------------------------------------------------------------------------
def basic_qc(adata: ad.AnnData, qc: dict) -> pd.DataFrame:
    """Per-nucleus QC metrics and filtering (the CELLxGENE release is already author-filtered;
    these thresholds only remove residual low-quality nuclei)."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    keep = (adata.obs["n_genes_by_counts"] >= qc["min_genes"]) & (adata.obs["pct_counts_mt"] <= qc["max_pct_mt"])
    summary = (
        adata.obs.assign(kept=keep)
        .groupby(["donor", "condition_label"], observed=True)
        .agg(n_nuclei=("kept", "size"), n_kept=("kept", "sum"),
             median_genes=("n_genes_by_counts", "median"), median_counts=("total_counts", "median"),
             median_pct_mt=("pct_counts_mt", "median"))
        .reset_index()
    )
    adata._inplace_subset_obs(keep.to_numpy())
    return summary


def composition(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cell-type proportions per donor and a Mann-Whitney test (test vs reference) per cell type."""
    counts = pd.crosstab(adata.obs["donor"], adata.obs["cell_type"])
    props = counts.div(counts.sum(axis=1), axis=0)
    donor_cond = adata.obs.groupby("donor", observed=True)["condition_label"].first()
    long = props.stack().rename("proportion").reset_index()
    long["condition"] = long["donor"].map(donor_cond)
    rows = []
    test_label = adata.obs.loc[adata.obs["condition"] == "test", "condition_label"].iloc[0]
    for ct, sub in long.groupby("cell_type"):
        a = sub.loc[sub["condition"] == test_label, "proportion"]
        b = sub.loc[sub["condition"] != test_label, "proportion"]
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) and len(b) else np.nan
        rows.append({"cell_type": ct, "median_test": a.median(), "median_reference": b.median(),
                     "n_test_donors": len(a), "n_reference_donors": len(b), "pvalue": p})
    stats = pd.DataFrame(rows)
    stats["fdr_bh"] = multipletests(stats["pvalue"].fillna(1), method="fdr_bh")[1]
    return long, stats


def aggregate_pseudobulk(adata: ad.AnnData, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts per donor x cell type; keep samples with at least `min_cells` nuclei."""
    groups = adata.obs.groupby(["cell_type", "donor"], observed=True).indices
    X = adata.X.tocsr()
    rows, meta = [], []
    for (ct, donor), idx in groups.items():
        if len(idx) < min_cells:
            continue
        rows.append(np.asarray(X[idx].sum(axis=0)).ravel())
        cond = adata.obs["condition"].iloc[idx[0]]
        meta.append({"sample": f"{ct}|{donor}", "cell_type": ct, "donor": donor,
                     "condition": cond, "condition_label": adata.obs["condition_label"].iloc[idx[0]],
                     "n_nuclei": len(idx)})
    if not meta:
        raise ValueError(f"No donor x cell-type sample has >= {min_cells} nuclei; lower pseudobulk.min_cells_per_sample")
    meta_df = pd.DataFrame(meta).set_index("sample")
    counts = pd.DataFrame(np.rint(np.vstack(rows)).astype(np.int64), index=meta_df.index, columns=adata.var_names)
    return counts, meta_df


def integrate_subsample(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """Donor-balanced subsample -> normalise, HVGs, PCA, Harmony on donor, UMAP (visualisation only)."""
    qc = cfg["qc"]
    rng = np.random.default_rng(int(cfg.get("seed", 42)))
    idx = []
    for _, sub_idx in adata.obs.groupby("donor", observed=True).indices.items():
        n = min(len(sub_idx), int(qc["nuclei_per_donor_for_umap"]))
        idx.extend(rng.choice(sub_idx, size=n, replace=False))
    sub = adata[np.sort(idx)].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.highly_variable_genes(sub, n_top_genes=int(qc["n_top_genes"]), batch_key="donor")
    sub = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=int(qc["n_pcs"]), random_state=int(cfg.get("seed", 42)))
    import harmonypy

    ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, "donor", random_state=int(cfg.get("seed", 42)), verbose=False)
    z = np.asarray(ho.Z_corr)
    sub.obsm["X_pca_harmony"] = z.T if z.shape[0] != sub.n_obs else z
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=int(qc["n_neighbors"]), random_state=int(cfg.get("seed", 42)))
    sc.tl.umap(sub, random_state=int(cfg.get("seed", 42)))
    return sub


def plot_fig1(sub: ad.AnnData, long: pd.DataFrame, stats: pd.DataFrame, out_path: str) -> None:
    fig = plt.figure(figsize=(16, 5.2), dpi=180)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25])
    um = sub.obsm["X_umap"]
    ax = fig.add_subplot(gs[0])
    cts = sorted(sub.obs["cell_type"].unique())
    cmap = plt.get_cmap("tab10")
    for i, ct in enumerate(cts):
        m = (sub.obs["cell_type"] == ct).to_numpy()
        ax.scatter(um[m, 0], um[m, 1], s=2, color=cmap(i % 10), label=ct, rasterized=True)
    ax.set_title("Harmony-integrated UMAP by published cell type", fontweight="bold", fontsize=10)
    ax.legend(markerscale=4, fontsize=6.5, frameon=False, loc="best")
    ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(gs[1])
    for lab, col in zip(sorted(sub.obs["condition_label"].unique()), ["#c05621", "#2b6cb0"]):
        m = (sub.obs["condition_label"] == lab).to_numpy()
        ax.scatter(um[m, 0], um[m, 1], s=2, color=col, alpha=0.6, label=lab, rasterized=True)
    ax.set_title("UMAP by cycle preparation", fontweight="bold", fontsize=10)
    ax.legend(markerscale=4, fontsize=8, frameon=False)
    ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(gs[2])
    order = stats.sort_values("median_reference", ascending=False)["cell_type"].tolist()
    labels = sorted(long["condition"].unique())
    width = 0.38
    for j, (lab, col) in enumerate(zip(labels, ["#c05621", "#2b6cb0"])):
        for i, ct in enumerate(order):
            vals = long[(long["cell_type"] == ct) & (long["condition"] == lab)]["proportion"]
            ax.boxplot(vals, positions=[i + (j - 0.5) * width], widths=width * 0.9, patch_artist=True,
                       boxprops=dict(facecolor=col, alpha=0.55), medianprops=dict(color="black"),
                       showfliers=False)
            ax.scatter(np.full(len(vals), i + (j - 0.5) * width), vals, s=8, color=col, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([o.replace(" of endometrium", "").replace(" cell", "") for o in order], rotation=40, ha="right", fontsize=7)
    ax.set_ylabel("Fraction of nuclei per donor")
    ax.set_title("Cell-type composition per donor", fontweight="bold", fontsize=10)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, alpha=0.55, label=l) for l, c in zip(labels, ["#c05621", "#2b6cb0"])], fontsize=8, frameon=False)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_step(cfg: dict, out: dict) -> None:
    d = cfg["dataset"]
    path = download_dataset(d["url"], Path(d["h5ad_path"]))
    adata = load_counts(path, cfg)
    qc_summary = basic_qc(adata, cfg["qc"])
    long, comp_stats = composition(adata)
    counts, meta = aggregate_pseudobulk(adata, int(cfg["pseudobulk"]["min_cells_per_sample"]))
    sub = integrate_subsample(adata, cfg)

    for key, df in (("qc", qc_summary), ("props", long), ("comp_stats", comp_stats)):
        Path(out[key]).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out[key], index=False)
    counts.to_csv(out["pb_counts"], compression="gzip")
    meta.to_csv(out["pb_meta"])
    Path(out["manifest"]).write_text(json.dumps({
        "dataset_url": d["url"], "citation": d["citation"], "license": d["license"],
        "sha256": sha256(path), "n_nuclei_after_qc": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "n_donors": int(adata.obs["donor"].nunique()), "n_pseudobulk_samples": int(len(meta)),
    }, indent=2))
    plot_fig1(sub, long, comp_stats, out["fig1"])


if __name__ == "__main__":
    if "snakemake" in globals():
        run_step(snakemake.params.cfg, {k: getattr(snakemake.output, k) for k in snakemake.output.keys()})
