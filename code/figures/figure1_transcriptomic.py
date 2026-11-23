"""Figure 1. Transcriptomic derivation of the four microenvironmental inputs.
(a) selection workflow of the four genes from the TCGA-UCEC transcriptome; (b) multivariable Cox model of the four genes in
the training cohort; (c) Hallmark scores of the proliferative and inflammatory programmes in high-risk and low-risk
tumours; (d, e) overall survival by risk group in the training and testing cohorts; (f) share of high-risk tumours
inside and outside the seven input contexts that contain a proliferative attractor."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style as S
from paths import FIGURES, TAB

S.apply()
HIGH, LOW = S.RUST, S.BLUE


def panel_funnel(ax, W, H):
    """Selection workflow drawn with vertical connectors: counts in the boxes, the selection step beside each arrow."""
    import matplotlib.patches as mp
    from orthonet import draw_edge

    F = pd.read_csv(TAB + "/L1_selection_funnel.csv")
    G = pd.read_csv(TAB + "/G1_signature_coefficients.csv")
    boxes = ["519 TCGA tumours"] + [f"{v} genes" for v in F.genes[:-1]] + [", ".join(G.gene)]
    steps = ["differential expression", "immune-related", "univariate Cox", "penalised Cox", "stepwise Cox"]
    bw, bh, x = 1.30, 0.15, 0.70
    pitch = (H - bh - 0.02) / (len(boxes) - 1)
    for i, t in enumerate(boxes):
        y = H - 0.01 - bh / 2 - i * pitch
        last = i == len(boxes) - 1
        ax.add_patch(mp.FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh, boxstyle="round,pad=0,rounding_size=0.02",
                                       fc=S.INK if last else "white", ec=S.INK, lw=0.6))
        ax.text(x, y - 0.004, t, ha="center", va="center", fontsize=6.6 if last else 7.2,
                color="white" if last else S.INK)
        if i:
            draw_edge(ax, [(x, y + pitch - bh / 2), (x, y + bh / 2)], False, S.INK, lw=0.6, head=0.028)
            ax.text(x + 0.06, y + pitch / 2, steps[i - 1], ha="left", va="center", fontsize=6.8, color=S.GREY)


def panel_forest(ax):
    G = pd.read_csv(TAB + "/G1_signature_coefficients.csv")
    for i, r in G.iterrows():
        col = HIGH if r.HR > 1 else LOW
        ax.plot([r.HR_lower_95, r.HR_upper_95], [i, i], color=col, lw=1.1)
        ax.plot(r.HR, i, marker="s", ms=5, color=col)
        if r.HR > 1:
            ax.text(r.HR_upper_95 * 1.08, i, f"{r.HR:.2f}", va="center", fontsize=7.5)
        else:
            ax.text(r.HR_lower_95 / 1.08, i, f"{r.HR:.2f}", va="center", ha="right", fontsize=7.5)
    ax.axvline(1, color=S.GREY, lw=0.5, ls=(0, (3, 2)))
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.25", "0.5", "1", "2"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(0.15, 3.0)
    ax.set_yticks(range(len(G)))
    ax.set_yticklabels(G.gene, fontsize=8)
    ax.set_ylim(len(G) - 0.5, -0.5)
    ax.set_xlabel("Hazard ratio (95% CI)")


def panel_km(ax, cohort, title):
    C = pd.read_csv(TAB + "/L2_km_curves.csv").query("cohort == @cohort")
    T = pd.read_csv(TAB + "/L3_risk_group_tests.csv").set_index("cohort").loc[cohort]
    for g, col in (("low risk", LOW), ("high risk", HIGH)):
        d = C[C.group == g]
        ax.step(d.time, d.survival, where="post", color=col, lw=1.2, label=f"{g} (n = {int(d.n.iloc[0])})")
        ax.fill_between(d.time, d.lower, d.upper, step="post", color=col, alpha=0.10, lw=0)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Years")
    if cohort == "training":
        ax.set_ylabel("Overall survival")
    ax.set_title(title, fontsize=8.5, pad=3)
    p = T.log_rank_P
    ptxt = f"log-rank P = {S.sci(p)}" if p < 1e-3 else f"log-rank P = {p:.4f}"
    ax.text(0.04, 0.05, f"HR {T.HR_high_vs_low:.2f} ({T.HR_lower_95:.2f} to {T.HR_upper_95:.2f})\n" + ptxt,
            transform=ax.transAxes, fontsize=7.5, va="bottom")
    if cohort == "training":
        ax.legend(loc="upper right", fontsize=7, handlelength=1.4, labelspacing=0.2, frameon=True, framealpha=1, edgecolor="none")
    else:
        ax.text(0.97, 0.76, "low risk n = %d\nhigh risk n = %d" % (int(C[C.group == "low risk"].n.iloc[0]), int(C[C.group == "high risk"].n.iloc[0])),
                transform=ax.transAxes, ha="right", va="top", fontsize=7)


def panel_hallmark(ax):
    H = pd.read_csv(TAB + "/L4b_hallmark_scores_with_groups.csv")
    L = pd.read_csv(TAB + "/L4_hallmark_by_risk.csv").set_index("hallmark")
    sets = [("E2F_TARGETS", "E2F\ntargets"), ("MYC_TARGETS_V1", "MYC\ntargets"),
            ("IL6_JAK_STAT3_SIGNALING", "IL6-JAK-\nSTAT3"), ("EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT")]
    rng = np.random.default_rng(1)
    for i, (h, lab) in enumerate(sets):
        for g, col, off in (("Low", LOW, -0.18), ("High", HIGH, 0.18)):
            v = H.loc[H.Risk_group == g, h].dropna().values
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.scatter(i + off + rng.uniform(-0.07, 0.07, len(v)), v, s=2.0, color=col, alpha=0.25, lw=0)
            ax.plot([i + off, i + off], [q1, q3], color=col, lw=2.6, solid_capstyle="butt")
            ax.plot([i + off - 0.1, i + off + 0.1], [med, med], color=S.INK, lw=1.0)
        q = L.loc[h, "FDR"]
        txt = "<0.001" if q < 0.001 else f"{q:.3f}" if q < 0.01 else f"{q:.2f}"
        ax.text(i, 1.12, txt, ha="center", fontsize=6.8, color=S.INK if q < 0.05 else S.GREY)
    ax.text(-0.75, 1.12, "FDR", ha="right", fontsize=6.8, color=S.GREY)
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([s[1] for s in sets], fontsize=7.6)
    ax.set_xlim(-0.6, len(sets) - 0.4)
    ax.set_ylim(-1.15, 1.32)
    ax.set_ylabel("Hallmark score")
    ax.plot([], [], color=LOW, lw=2.6, label="low risk")
    ax.plot([], [], color=HIGH, lw=2.6, label="high risk")
    ax.legend(loc="lower left", fontsize=7, ncol=2, handlelength=1.2, columnspacing=0.8, bbox_to_anchor=(0, -0.04))


def panel_contexts(ax):
    R = pd.read_csv(TAB + "/L5_risk_by_proliferative_context.csv").set_index("proliferative_context")
    rows = [(True, "Proliferative\ncontexts"), (False, "Other\ncontexts")]
    for i, (k, lab) in enumerate(rows):
        hi, lo = R.loc[k, "High"], R.loc[k, "Low"]
        n = hi + lo
        ax.barh(i, hi / n, color=HIGH, height=0.55)
        ax.barh(i, lo / n, left=hi / n, color=LOW, height=0.55, alpha=0.35)
        ax.text(hi / n - 0.02, i, f"{hi}/{n}", ha="right", va="center", fontsize=7.2, color="white")
    ax.set_yticks([0, 1])
    ax.set_yticklabels([r[1] for r in rows], fontsize=7.6)
    ax.set_ylim(1.6, -0.6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of tumours at high risk")
    ax.set_xticks([0, 0.5, 1])
    p = R.fisher_P.iloc[0]
    ax.text(1.0, 1.06, f"Fisher P = {S.sci(p)}",
            transform=ax.transAxes, ha="right", fontsize=7.5)


if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 5.2))
    top, bot, h = 0.58, 0.08, 0.34
    W, H = 2.3, h * 5.2
    ax = fig.add_axes([0.0, top, W / 7.2, h]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    panel_funnel(ax, W, H); S.panel(ax, "a", dx=0.04, dy=1.03)
    ax = fig.add_axes([0.47, top, 0.17, h]); panel_forest(ax); S.panel(ax, "b", dx=-0.5, dy=1.03)
    ax = fig.add_axes([0.745, top, 0.25, h]); panel_hallmark(ax); S.panel(ax, "c", dx=-0.2, dy=1.03)
    ax = fig.add_axes([0.075, bot, 0.255, h]); panel_km(ax, "training", "Training cohort"); S.panel(ax, "d", dx=-0.2, dy=1.03)
    ax = fig.add_axes([0.405, bot, 0.255, h]); panel_km(ax, "testing", "Testing cohort"); S.panel(ax, "e", dx=-0.2, dy=1.03)
    ax = fig.add_axes([0.83, bot + 0.07, 0.165, 0.2]); panel_contexts(ax); S.panel(ax, "f", dx=-0.7, dy=1.12)
    S.save(fig, FIGURES, "Figure1")
