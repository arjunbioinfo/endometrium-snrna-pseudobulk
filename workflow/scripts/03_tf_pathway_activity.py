"""Step 3: transcription-factor and signalling-pathway activity behind the pseudobulk DE results.

For every tested cell type, the PyDESeq2 Wald statistics (programmed vs natural cycle) are scored
against prior-knowledge networks with decoupler:

  * TF activity: CollecTRI regulons (signed TF -> target links), univariate linear model (ULM).
  * Pathway activity: PROGENy footprints (top responsive genes per pathway), multivariate linear
    model (MLM).

Positive scores mean the regulator/pathway targets shift towards higher expression in the test
condition. Hormone receptors (PGR, ESR1, AR, NR3C1) are reported separately because the two cycle
preparations differ in steroid-hormone exposure.

References: Badia-i-Mompel et al., Bioinformatics Advances 2022 (decoupler); Müller-Dott et al.,
Nucleic Acids Research 2023 (CollecTRI); Schubert et al., Nature Communications 2018 (PROGENy).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

HORMONE_RECEPTORS = ["PGR", "ESR1", "AR", "NR3C1"]


def load_network(kind: str, cache: Path, top: int = 500) -> pd.DataFrame:
    """CollecTRI or PROGENy from OmniPath, cached as CSV so later runs work offline."""
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.stat().st_size > 0:
        return pd.read_csv(cache)
    import decoupler as dc

    if kind == "collectri":
        net = dc.get_collectri(organism="human", split_complexes=False)[["source", "target", "weight"]]
    elif kind == "progeny":
        net = dc.get_progeny(organism="human", top=top)[["source", "target", "weight"]]
    else:
        raise ValueError(kind)
    net.to_csv(cache, index=False)
    return net


def score_contrasts(de: pd.DataFrame, net: pd.DataFrame, method: str, min_targets: int) -> pd.DataFrame:
    """Score each cell type's Wald statistics against `net`; returns long table with BH FDR per cell type."""
    import decoupler as dc

    frames = []
    for ct, sub in de.dropna(subset=["stat"]).groupby("cell_type"):
        mat = sub.drop_duplicates("gene").set_index("gene")["stat"].to_frame(ct).T
        runner = dc.run_ulm if method == "ulm" else dc.run_mlm
        est, pvals = runner(mat=mat, net=net, source="source", target="target", weight="weight",
                            min_n=min_targets, verbose=False)
        df = pd.DataFrame({"source": est.columns, "score": est.iloc[0].to_numpy(), "pvalue": pvals.iloc[0].to_numpy()})
        df.insert(0, "cell_type", ct)
        df["fdr_bh"] = multipletests(df["pvalue"].fillna(1), method="fdr_bh")[1]
        df["n_targets_tested"] = df["source"].map(
            net[net["target"].isin(mat.columns)].groupby("source").size()).fillna(0).astype(int)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _short(ct: str) -> str:
    return ct.replace(" of endometrium", "").replace(" of uterus", "").replace(" cell", "")


def _heatmap(ax, piv: pd.DataFrame, sig: pd.DataFrame, title: str, cbar_label: str, fig) -> None:
    vmax = float(np.nanmax(np.abs(piv.to_numpy()))) if piv.size else 1.0
    im = ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([_short(c) for c in piv.columns], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=8)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if bool(sig.iat[i, j]):
                ax.text(j, i, "*", ha="center", va="center", fontsize=9, color="black")
    ax.set_title(title, fontweight="bold", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbar_label)


def plot_fig3(tf: pd.DataFrame, pw: pd.DataFrame, cfg: dict, out_path: str) -> None:
    d = cfg["dataset"]
    n_top = int(cfg["activity"]["n_top_tfs"])
    fig, axes = plt.subplots(1, 3, figsize=(17, 7.2), dpi=180, gridspec_kw={"width_ratios": [1.15, 0.9, 1.0]})
    if not tf.empty:
        top = (tf[tf["fdr_bh"] < 0.05].assign(a=lambda x: x["score"].abs())
               .sort_values("a", ascending=False)["source"].drop_duplicates().head(n_top).tolist())
        if top:
            piv = tf[tf["source"].isin(top)].pivot(index="source", columns="cell_type", values="score").reindex(top)
            sig = tf[tf["source"].isin(top)].pivot(index="source", columns="cell_type", values="fdr_bh").reindex(top) < 0.05
            _heatmap(axes[0], piv, sig, f"Top transcription-factor activities (CollecTRI)", "ULM score", fig)
        hr = tf[tf["source"].isin(HORMONE_RECEPTORS)]
        piv = hr.pivot(index="source", columns="cell_type", values="score").reindex([h for h in HORMONE_RECEPTORS if h in set(hr["source"])])
        sig = hr.pivot(index="source", columns="cell_type", values="fdr_bh").reindex(piv.index) < 0.05
        _heatmap(axes[1], piv, sig, "Steroid-hormone receptors", "ULM score", fig)
    if not pw.empty:
        piv = pw.pivot(index="source", columns="cell_type", values="score")
        sig = pw.pivot(index="source", columns="cell_type", values="fdr_bh").reindex(index=piv.index, columns=piv.columns) < 0.05
        _heatmap(axes[2], piv, sig, "Pathway activities (PROGENy)", "MLM score", fig)
    fig.suptitle(f"Positive = targets shifted towards {d['test_level']} (programmed) vs {d['reference_level']} cycle;  * FDR < 0.05",
                 fontsize=9.5)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_step(de_path: str, cfg: dict, out_tf: str, out_pw: str, out_fig: str) -> None:
    a = cfg["activity"]
    ref = Path(a.get("reference_dir", "data/reference"))
    de = pd.read_csv(de_path)
    collectri = load_network("collectri", ref / "collectri_human.csv")
    progeny = load_network("progeny", ref / f"progeny_human_top{int(a['progeny_top'])}.csv", top=int(a["progeny_top"]))
    tf = score_contrasts(de, collectri, "ulm", int(a["min_targets"]))
    pw = score_contrasts(de, progeny, "mlm", int(a["min_targets"]))
    for df, p in ((tf, out_tf), (pw, out_pw)):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)
    plot_fig3(tf, pw, cfg, out_fig)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_step(snakemake.input.de, snakemake.params.cfg, snakemake.output.tf, snakemake.output.pathways,
                 snakemake.output.fig3)
