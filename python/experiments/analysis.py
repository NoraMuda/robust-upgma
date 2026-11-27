"""Every table, figure and quoted number of the manuscript, from the trial-level
results in results/python/*.csv.gz.

    python analysis.py
Output: results/python/tables/*.csv, results/python/figures/*.png,
        results/python/key_numbers.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from common import FCPS, REAL, RESULTS as RES

OUT = RES
TAB, FIG = OUT / "tables", OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

ALL = ["R-UPGMA", "W-UPGMA", "UPGMA", "UPGMA+cut", "Complete", "Complete+cut", "Single", "PAM", "HDBSCAN"]
MAIN = ["R-UPGMA", "W-UPGMA", "UPGMA", "UPGMA+cut", "Complete", "Complete+cut", "PAM", "HDBSCAN"]
FIGM = ["R-UPGMA", "UPGMA", "UPGMA+cut", "Complete+cut", "PAM", "HDBSCAN"]
COLORS = dict(zip(FIGM, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]))
SEP = ["well", "moderate", "overlap"]
SEP_LAB = {"well": "Well separated (Δ=6)", "moderate": "Moderate (Δ=4)", "overlap": "Overlapping (Δ=2.5)"}
CORR = ["none", "salt_pepper", "amplification", "mixed", "gaussian", "structured", "row"]
CORR_LAB = {"none": "None", "salt_pepper": "Salt-and-pepper", "amplification": "Amplification",
            "mixed": "Mixed", "gaussian": "Gaussian", "structured": "Cluster-structured",
            "row": "Observation outliers"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#888", "axes.linewidth": 0.6,
                     "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.5, "axes.axisbelow": True})
KEY: dict = {}


def q(s, a):
    return float(np.quantile(s, a))


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = run
    return adj


def load_sim():
    df = pd.read_csv(RES / "sim_trials.csv.gz")
    return df[df.method.isin(ALL)]


def simulation():
    df = load_sim()
    KEY["sim_n"] = int(df[["cond", "rep"]].drop_duplicates().shape[0])
    s1 = (df.groupby(["structure", "separation", "corruption", "level", "method"]).ari
          .agg(n="size", mean="mean", sd="std", q05=lambda s: q(s, .05), q10=lambda s: q(s, .1),
               median="median", q90=lambda s: q(s, .9), q95=lambda s: q(s, .95),
               fail_rate=lambda s: np.mean(s < .5)).reset_index())
    s1.round(4).to_csv(TAB / "S1_sim_ari_by_condition_and_method.csv", index=False)

    ov = df.groupby("method").ari.agg(mean="mean", median="median", fail=lambda s: np.mean(s < .5))
    ov.loc[ALL].round(4).to_csv(TAB / "sim_overall.csv")
    KEY["sim_overall"] = ov.loc[ALL].round(3).to_dict(orient="index")

    sub = df[df.level.isin([0.0, 0.10])]
    t4 = sub.groupby(["separation", "corruption", "method"]).ari.median().unstack()[MAIN].reset_index()
    t4["separation"] = pd.Categorical(t4.separation, SEP)
    t4["corruption"] = pd.Categorical(t4.corruption, CORR)
    t4 = t4.sort_values(["separation", "corruption"])
    t4.round(2).to_csv(TAB / "T4_median_ari_eps10.csv", index=False)
    fr = sub.groupby(["separation", "corruption", "method"]).ari.apply(lambda s: np.mean(s < .5)).unstack()[MAIN]
    fr.round(3).to_csv(TAB / "T4b_failure_eps10.csv")
    KEY["median_eps10"] = {f"{a}|{b}|{c}": round(float(v), 2) for (a, b, c), v in
                           sub.groupby(["separation", "corruption", "method"]).ari.median().items()}

    p = df.pivot_table(index=["structure", "separation", "corruption", "level", "rep"], columns="method",
                       values="ari").reset_index()
    rows = []
    for comp in ["PAM", "UPGMA+cut"]:
        for (sep, cor), g in p.groupby(["separation", "corruption"]):
            d = (g["R-UPGMA"] - g[comp]).to_numpy()
            try:
                pv = wilcoxon(d, zero_method="zsplit").pvalue
            except ValueError:
                pv = 1.0
            rows.append(dict(comparator=comp, separation=sep, corruption=cor, n=len(d), mean=d.mean(),
                             median=np.median(d), better=np.mean(d > .05), worse=np.mean(d < -.05), p=pv))
    t5 = pd.DataFrame(rows)
    t5["p_holm"] = np.nan
    for comp in ["PAM", "UPGMA+cut"]:
        m = t5.comparator == comp
        t5.loc[m, "p_holm"] = holm(t5.loc[m, "p"].values)
    t5["separation"] = pd.Categorical(t5.separation, SEP)
    t5["corruption"] = pd.Categorical(t5.corruption, CORR)
    t5 = t5.sort_values(["comparator", "separation", "corruption"])
    t5.round(4).to_csv(TAB / "T5_paired_vs_PAM_and_UPGMAcut.csv", index=False)
    for comp in ["PAM", "UPGMA+cut", "UPGMA", "W-UPGMA", "Complete+cut"]:
        d = p["R-UPGMA"] - p[comp]
        KEY[f"paired_vs_{comp}"] = dict(mean=round(float(d.mean()), 3), better=round(float(np.mean(d > .05)), 3),
                                        worse=round(float(np.mean(d < -.05)), 3))
    dc = (p["R-UPGMA"] - p["UPGMA"])[p.corruption == "none"]
    KEY["clean_vs_UPGMA"] = dict(mean=round(float(dc.mean()), 3), better=round(float(np.mean(dc > .05)), 3),
                                 worse=round(float(np.mean(dc < -.05)), 3))

    # Figure 1
    fig, axes = plt.subplots(len(CORR), 3, figsize=(7.2, 9.6), sharey=True)
    for i, cor in enumerate(CORR):
        for j, sep in enumerate(SEP):
            ax = axes[i, j]
            g = sub[(sub.corruption == cor) & (sub.separation == sep)]
            data = [g[g.method == m].ari.to_numpy() for m in FIGM]
            parts = ax.violinplot(data, showextrema=False, widths=0.8, bw_method=0.15)
            for b, mth in zip(parts["bodies"], FIGM):
                b.set_facecolor(COLORS[mth]); b.set_edgecolor("white"); b.set_linewidth(0.6); b.set_alpha(0.85)
            ax.scatter(range(1, 7), [np.median(x) for x in data], s=9, color="#222", zorder=3, lw=0)
            ax.set_ylim(-0.1, 1.05)
            ax.set_xticks(range(1, 7))
            ax.set_xticklabels(FIGM if i == len(CORR) - 1 else [], rotation=55, ha="right", fontsize=7)
            if i == 0:
                ax.set_title(SEP_LAB[sep], fontsize=8)
            if j == 0:
                ax.set_ylabel(CORR_LAB[cor] + ("" if cor == "none" else " 10%"), fontsize=7)
    fig.supylabel("Adjusted Rand Index (per simulated dataset)", fontsize=8)
    fig.tight_layout(h_pad=0.4, w_pad=0.4)
    fig.savefig(FIG / "Fig1_sim_ari_distributions.png", dpi=300)
    plt.close(fig)

    # Figure 2: ECDF of paired differences vs PAM
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
    seq = ["#0d366b", "#256abf", "#5598e7"]
    styles = {"none": ("#999999", "--"), "salt_pepper": (seq[0], "-"), "amplification": (seq[1], "-"),
              "mixed": (seq[2], "-"), "structured": ("#eb6834", "-"), "row": ("#444444", ":")}
    for j, sep in enumerate(SEP):
        ax = axes[j]
        for cor, (c, ls) in styles.items():
            g = p[(p.separation == sep) & (p.corruption == cor)]
            d = np.sort((g["R-UPGMA"] - g["PAM"]).to_numpy())
            ax.step(d, np.arange(1, len(d) + 1) / len(d), where="post", color=c, ls=ls, lw=1.4, label=CORR_LAB[cor])
        ax.axvline(0, color="#888", lw=0.6)
        ax.set_xlim(-0.8, 1.0)
        ax.set_title(SEP_LAB[sep], fontsize=8)
        ax.set_xlabel("ΔARI = ARI(R-UPGMA) − ARI(PAM)")
    axes[0].set_ylabel("Cumulative proportion of datasets")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, loc="lower center", ncol=6)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(FIG / "Fig2_sim_paired_vs_PAM_ecdf.png", dpi=300)
    plt.close(fig)


def ablation_and_dev():
    order = ["R-UPGMA", "w/o Step 1 (triangle)", "w/o Step 2 (MADn bounds)", "w/o Step 3 (trimming)",
             "w/o Step 4 (cut)"]
    df = pd.read_csv(RES / "sim_trials.csv.gz")
    df = df[(df.rep < 30) & df.method.isin(order)]     # ablations were run on replicates 0-29
    t = df.groupby(["corruption", "method"]).ari.mean().unstack()[order]
    t = t.loc[CORR]
    tot = df.groupby("method").ari.agg(mean="mean", fail=lambda s: np.mean(s < .5))[["mean", "fail"]].loc[order]
    t.loc["all (mean)"] = tot["mean"]
    t.loc["all (failure rate)"] = tot["fail"]
    t.round(3).to_csv(TAB / "T6_ablation.csv")
    KEY["ablation_overall"] = tot.round(3).to_dict(orient="index")
    KEY["ablation_by_corruption"] = t.round(3).to_dict()
    dev = pd.read_csv(RES / "development_round2.csv.gz")
    dv = dev.groupby("method").ari.agg(mean="mean", fail=lambda s: np.mean(s < .5)).round(3)
    dv.to_csv(TAB / "S5_development_variants.csv")
    KEY["development"] = dv.to_dict(orient="index")


def benchmarks():
    df = pd.read_csv(RES / "bench_trials.csv.gz")
    sets = FCPS + list(REAL)
    med = df.groupby(["dataset", "setting", "method"]).ari.median().unstack()[ALL]
    rows = []
    for ds in sets:
        for st in ["clean", "salt_pepper10", "mixed10", "row05"]:
            r = med.loc[(ds, st)]
            g = df[(df.dataset == ds) & (df.setting == st)]
            rows.append(dict(dataset=ds, setting=st, n=int(g.n.iloc[0]), K=int(g.K.iloc[0]), **r.round(2).to_dict()))
    t8 = pd.DataFrame(rows)
    t8.to_csv(TAB / "T8_benchmarks_median_ari.csv", index=False)
    full = df.groupby(["dataset", "setting", "method"]).agg(
        n_runs=("ari", "size"), ari_median=("ari", "median"), ari_q25=("ari", lambda s: q(s, .25)),
        ari_q75=("ari", lambda s: q(s, .75)), acc_median=("acc", "median"), sil_median=("sil", "median"),
        largest_share=("largest", "median")).round(4)
    full.to_csv(TAB / "S2_benchmarks_full.csv")
    df.groupby(["dataset", "setting", "method"]).sil.median().unstack()[ALL].round(3).to_csv(
        TAB / "S3_benchmarks_silhouette_clean_geometry.csv")
    m = med[MAIN]
    top = m.max(axis=1)
    KEY["bench_best_or_tied"] = {k: int((m[k] >= top - .01).sum()) for k in MAIN}
    KEY["bench_cells"] = int(len(m))
    dpam = m["R-UPGMA"] - m["PAM"]
    KEY["bench_vs_PAM"] = dict(min=round(float(dpam.min()), 2), max=round(float(dpam.max()), 2),
                               better=int((dpam > .05).sum()), worse=int((dpam < -.05).sum()))
    dr = m["R-UPGMA"] - m["UPGMA+cut"]
    KEY["bench_vs_UPGMAcut"] = dict(min=round(float(dr.min()), 2), max=round(float(dr.max()), 2),
                                     better=int((dr > .05).sum()), worse=int((dr < -.05).sum()))
    KEY["bench_median"] = {f"{d}|{s}|{k}": float(round(v, 2)) for (d, s), r in m.iterrows() for k, v in r.items()}

    nc = 5
    nr = int(np.ceil(len(sets) / nc))
    for st, fname in [("salt_pepper10", "Fig3_bench_saltpepper10.png"), ("mixed10", "Fig4_bench_mixed10.png")]:
        fig, axes = plt.subplots(nr, nc, figsize=(7.2, 1.95 * nr + 0.4), sharey=True)
        for ax, ds in zip(axes.ravel(), sets):
            g = df[(df.dataset == ds) & (df.setting == st)]
            data = [g[g.method == mm].ari.to_numpy() for mm in FIGM]
            parts = ax.violinplot(data, showextrema=False, widths=0.8, bw_method=0.3)
            for b, mth in zip(parts["bodies"], FIGM):
                b.set_facecolor(COLORS[mth]); b.set_edgecolor("white"); b.set_alpha(0.85)
            clean = med.loc[(ds, "clean")]
            ax.scatter(range(1, 7), [clean[mm] for mm in FIGM], marker="_", s=60, color="#222", lw=1.2, zorder=3)
            ax.scatter(range(1, 7), [np.median(x) for x in data], s=8, color="#222", zorder=3, lw=0)
            ax.set_title(f"{ds}\n(n={int(g.n.iloc[0])}, K={int(g.K.iloc[0])})", fontsize=7)
            ax.set_xticks(range(1, 7)); ax.set_ylim(-0.1, 1.05)
            ax.set_xticklabels(FIGM, rotation=60, ha="right", fontsize=6)
        for ax in axes.ravel()[len(sets):]:
            ax.axis("off")
        fig.supylabel("ARI (30 corrupted realisations)", fontsize=8)
        fig.tight_layout(h_pad=0.6)
        fig.savefig(FIG / fname, dpi=300)
        plt.close(fig)


def null():
    df = pd.read_csv(RES / "null_trials.csv.gz")
    t = df.groupby(["dataset", "setting"]).agg(
        dip_p=("dip_p", "median"), reject_raw=("dip_p", lambda s: np.mean(s < .05)),
        reject_wins=("dip_p_winsorised", lambda s: np.mean(s < .05)),
        largest=("largest", "median"), smallest=("smallest", "median"), sil=("sil", "median")).round(3)
    t.to_csv(TAB / "T7_null.csv")
    KEY["null"] = {f"{a}|{b}": r.to_dict() for (a, b), r in t.iterrows()}


def timing():
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist
    from robust_upgma import robust_upgma
    rng = np.random.default_rng(1)
    rows = []
    for n in [150, 1000, 2000, 4000]:
        d = pdist(rng.standard_normal((n, 10)))
        a, b = [], []
        for _ in range(3 if n < 4000 else 1):
            t0 = time.perf_counter(); linkage(d, "average"); a.append(time.perf_counter() - t0)
            t0 = time.perf_counter(); robust_upgma(d, 3); b.append(time.perf_counter() - t0)
        rows.append(dict(n=n, pairs=len(d), upgma_s=round(float(np.median(a)), 3),
                         rupgma_s=round(float(np.median(b)), 2)))
    pd.DataFrame(rows).to_csv(TAB / "T9_timing.csv", index=False)
    KEY["timing"] = rows


if __name__ == "__main__":
    simulation()
    ablation_and_dev()
    benchmarks()
    null()
    timing()
    (OUT / "key_numbers.json").write_text(json.dumps(KEY, indent=1, default=float))
    print(json.dumps({k: KEY[k] for k in KEY if not k.startswith(("median_eps10", "bench_median"))},
                     indent=1, default=float))
