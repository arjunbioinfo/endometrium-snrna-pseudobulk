"""Offline unit tests: no downloads, small synthetic inputs only."""

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
CFG = {
    "seed": 1,
    "dataset": {"symbol_column": "feature_name", "donor_column": "donor_id", "condition_column": "CycleType",
                "celltype_column": "cell_type", "test_level": "Controlled", "reference_level": "Natural"},
    "qc": {"min_genes": 1, "max_pct_mt": 50.0},
    "pseudobulk": {"min_cells_per_sample": 5, "min_donors_per_group": 3, "min_gene_count": 5,
                   "padj_threshold": 0.05, "lfc_threshold": 0.5},
}


def _load(rel_path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _toy_h5ad(tmp_path: Path, n_donors_per_group: int = 4, cells: int = 30, genes: int = 40) -> Path:
    """CELLxGENE-like h5ad: Ensembl IDs as var_names, symbols in raw.var['feature_name'],
    raw counts in .raw.X, normalised values in .X. Gene G0 is 4x higher in 'Controlled'."""
    rng = np.random.default_rng(0)
    obs, mats = [], []
    for g, cond in enumerate(["Controlled", "Natural"]):
        for d in range(n_donors_per_group):
            for ct in ["stromal", "epithelial"]:
                lam = np.full(genes, 5.0)
                if cond == "Controlled":
                    lam[0] = 20.0
                mats.append(rng.poisson(lam, size=(cells, genes)))
                obs += [{"donor_id": f"D{g}{d}", "CycleType": cond, "cell_type": ct}] * cells
    X = np.vstack(mats).astype(np.float32)
    var = pd.DataFrame({"feature_name": [f"G{i}" for i in range(genes - 1)] + ["MT-CO1"]},
                       index=[f"ENSG{i:011d}" for i in range(genes)])
    a = ad.AnnData(X=sp.csr_matrix(np.log1p(X)), obs=pd.DataFrame(obs, index=[f"c{i}" for i in range(len(obs))]), var=var)
    a.raw = ad.AnnData(X=sp.csr_matrix(X), obs=a.obs, var=var)
    p = tmp_path / "toy.h5ad"
    a.write_h5ad(p)
    return p


def test_load_counts_uses_raw_counts_and_symbols(tmp_path: Path):
    m1 = _load("workflow/scripts/01_download_qc_pseudobulk.py", "s1")
    a = m1.load_counts(_toy_h5ad(tmp_path), CFG)
    assert a.var_names[0] == "G0" and "MT-CO1" in a.var_names
    assert np.allclose(a.X.data, np.round(a.X.data))  # integer counts, not log-normalised
    assert set(a.obs["condition"]) == {"test", "reference"}
    assert set(a.obs["condition_label"]) == {"Controlled", "Natural"}


def test_pseudobulk_aggregation_sums_counts(tmp_path: Path):
    m1 = _load("workflow/scripts/01_download_qc_pseudobulk.py", "s1")
    a = m1.load_counts(_toy_h5ad(tmp_path), CFG)
    total = a.X.sum()
    counts, meta = m1.aggregate_pseudobulk(a, min_cells=5)
    assert len(meta) == 2 * 4 * 2 and (meta["n_nuclei"] == 30).all()
    assert counts.values.sum() == total
    import pytest
    with pytest.raises(ValueError, match="No donor x cell-type sample"):
        m1.aggregate_pseudobulk(a, min_cells=31)  # every sample too small


def test_pseudobulk_deseq2_detects_planted_gene(tmp_path: Path):
    m1 = _load("workflow/scripts/01_download_qc_pseudobulk.py", "s1")
    m2 = _load("workflow/scripts/02_pseudobulk_de.py", "s2")
    a = m1.load_counts(_toy_h5ad(tmp_path), CFG)
    counts, meta = m1.aggregate_pseudobulk(a, min_cells=5)
    de, summary = m2.run_all(counts, meta, CFG)
    top = de[de["cell_type"] == "stromal"].sort_values("padj").iloc[0]
    assert top["gene"] == "G0" and top["log2FoldChange"] > 1.5 and top["padj"] < 1e-6
    assert summary["tested"].all() and (summary["n_up"] >= 1).all()


def test_composition_stats_columns(tmp_path: Path):
    m1 = _load("workflow/scripts/01_download_qc_pseudobulk.py", "s1")
    a = m1.load_counts(_toy_h5ad(tmp_path), CFG)
    long, stats = m1.composition(a)
    assert set(stats["cell_type"]) == {"stromal", "epithelial"}
    assert {"pvalue", "fdr_bh", "median_test", "median_reference"} <= set(stats.columns)
    assert np.allclose(long.groupby("donor")["proportion"].sum(), 1.0)


# ---------------------------------------------------------------------------------------
# Step 3: TF and pathway activity (toy network, no download)
# ---------------------------------------------------------------------------------------
def test_tf_activity_sign_follows_targets():
    m3 = _load("workflow/scripts/03_tf_pathway_activity.py", "s3")
    genes = [f"G{i}" for i in range(60)]
    stat = np.zeros(60)
    stat[:10] = 6.0    # targets activated by TF_UP
    stat[10:20] = -6.0  # targets activated by TF_DOWN
    stat[20:] = np.random.default_rng(0).normal(0, 1, 40)
    de = pd.DataFrame({"cell_type": "stromal", "gene": genes, "stat": stat})
    net = pd.DataFrame({"source": ["TF_UP"] * 10 + ["TF_DOWN"] * 10,
                        "target": genes[:10] + genes[10:20], "weight": 1.0})
    res = m3.score_contrasts(de, net, "ulm", min_targets=5).set_index("source")
    assert res.loc["TF_UP", "score"] > 0 and res.loc["TF_DOWN", "score"] < 0
    assert res.loc["TF_UP", "fdr_bh"] < 0.05 and res.loc["TF_UP", "n_targets_tested"] == 10


def test_network_cache_is_reused(tmp_path: Path):
    m3 = _load("workflow/scripts/03_tf_pathway_activity.py", "s3")
    cache = tmp_path / "collectri_human.csv"
    pd.DataFrame({"source": ["A"], "target": ["B"], "weight": [1]}).to_csv(cache, index=False)
    net = m3.load_network("collectri", cache)  # must not download when the cache exists
    assert net.to_dict("list") == {"source": ["A"], "target": ["B"], "weight": [1]}


# ---------------------------------------------------------------------------------------
# Step 4: leave-one-donor-out robustness
# ---------------------------------------------------------------------------------------
def test_lodo_retains_planted_gene(tmp_path: Path):
    m1 = _load("workflow/scripts/01_download_qc_pseudobulk.py", "s1")
    m2 = _load("workflow/scripts/02_pseudobulk_de.py", "s2")
    m4 = _load("workflow/scripts/04_leave_one_donor_out.py", "s4")
    a = m1.load_counts(_toy_h5ad(tmp_path), CFG)
    counts, meta = m1.aggregate_pseudobulk(a, min_cells=5)
    meta = meta[meta["cell_type"] == "stromal"]
    de, _ = m2.run_all(counts.loc[meta.index], meta, CFG)
    cfg = {**CFG, "robustness": {"min_degs": 1}}
    genes, summary = m4.leave_one_donor_out(counts.loc[meta.index], meta, de, cfg)
    g0 = genes.set_index("gene").loc["G0"]
    assert g0["n_runs"] == 8 and g0["fraction_retained"] == 1.0
    assert summary.iloc[0]["retained_all_runs"] >= 1
