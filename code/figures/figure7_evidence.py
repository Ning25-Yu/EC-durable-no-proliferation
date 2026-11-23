"""Figure 7. Molecular and functional evidence from endometrial cells, endometrial tumours and published withdrawal
experiments.
(a) Ishikawa cells (GSE236228), control or penfluridol: Proliferation and Apoptosis margins of each sample;
(b) AN3-CA cells (GSE61481), aFGF or aFGF plus AZD4547: Proliferation margins; (c) CPTAC-UCEC: p27 protein in tumours
of the seven proliferative contexts and of the other contexts; rows 2 and 3 as follows, lettered (d) to (i).
(d) model: durable arrest probability of a p16 pulse against exposure duration; (e, f) Daigh et al. 2024, MCF10A:
phospho-Rb and EdU positivity 24 h and 72 h after removal of p16 that had been induced for 1, 3 or 5 days;
(g) Daigh et al.: re-entry after release with p21, p27 or p53 depleted, relative to non-targeting siRNA;
(h) Paulsohn et al. 2026, MIA PaCa-2: growth after withdrawal of palbociclib with sotorasib (day 3), control and
CDKN1B-depleted cells, normalised to day 3; (i) post-release fold-change ratios of p27-depleted or knockout cells
relative to controls. The murine readouts are cell confluence and cell viability, not tumour measurements.
The human-cell endpoint is day 10 and the murine-cell endpoint is day 7."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style as S
from paths import FIGURES, P16_TAB, PANC_TAB, TAB

S.apply()
SRC = pd.read_csv(P16_TAB + "/all_pulse_source_values.csv")


def panel_model(ax):
    K = pd.read_csv(TAB + "/K2_duration_response.csv").query("pulse == 'CDKN2A ON'")
    ax.plot(K.L, K.arrest, color=S.BLUE, lw=1.3)
    for L in (1, 3, 5):
        y = float(np.interp(L, K.L, K.arrest))
        ax.plot(L, y, marker="o", ms=4.5, mfc="white", mec=S.BLUE, mew=1.0, zorder=3)
    ax.set_xlim(0, 8)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("p16 exposure (model time)")
    ax.set_ylabel("Durable arrest probability")


def panel_duration(ax, readout, title, ymax):
    d = SRC[(SRC.readout == readout) & (SRC.role == "pulse_duration")]
    cyc = SRC[(SRC.readout == readout) & (SRC.role == "cycling_control")]
    for hours, col, mfc, off in ((24, S.BLUE, S.BLUE, -0.08), (72, S.RUST, "white", 0.08)):
        g = d[d.nominal_readout_hours == hours]
        m = g.groupby("pulse_duration_days").positive_cells_percent.agg(["mean", "std"])
        ax.errorbar(m.index + off, m["mean"], yerr=m["std"], color=col, lw=1.1, marker="o", ms=4.5, mfc=mfc, mec=col, mew=1.0,
                    capsize=2, label=f"{hours} h after removal", zorder=3)
        ax.scatter(g.pulse_duration_days + off, g.positive_cells_percent, s=6, color=col, alpha=0.5, lw=0, zorder=2)
        c = cyc[cyc.nominal_readout_hours == hours].positive_cells_percent.mean()
        ax.axhline(c, color=col, lw=0.5, ls=(0, (3, 2)), zorder=1)
    ax.text(5.35, cyc.positive_cells_percent.max() + 0.02 * ymax, "cycling controls", ha="right", va="bottom", fontsize=6.8, color=S.GREY)
    ax.set_xticks([1, 3, 5])
    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("p16 induction (days)")
    ax.set_ylabel(title)


def panel_depletion(ax):
    D = pd.read_csv(P16_TAB + "/depletion_effect_estimates.csv").set_index("condition")
    R = pd.read_csv(P16_TAB + "/depletion_source_ratios.csv")
    order = ["si-p21", "si-p53", "si-p27"]
    for i, c in enumerate(order):
        r = D.loc[c]
        col = S.INK if r.log_t95_lower > 1 else S.GREY
        ax.plot([r.log_t95_lower, r.log_t95_upper], [i, i], color=col, lw=1.0)
        ax.plot(r.geometric_mean_ratio, i, marker="o", ms=5, mfc=col, mec=col, zorder=3)
        v = R[R.condition == c].reentry_ratio_to_Nts
        ax.scatter(v, [i + 0.22] * len(v), s=7, color=col, alpha=0.6, lw=0)
    ax.axvline(1, color=S.GREY, lw=0.5, ls=(0, (3, 2)))
    ax.set_yticks(range(3))
    ax.set_yticklabels(["p21 depleted", "p53 depleted", "p27 depleted"], fontsize=8)
    ax.set_ylim(2.6, -0.5)
    ax.set_xlim(0.9, 1.3)
    ax.set_xlabel("Re-entry relative to control siRNA")


def panel_traj(ax):
    T = pd.read_csv(PANC_TAB + "/all_trajectories_and_normalizations.csv")
    T = T[(T.panel == "G") & (T.drug == "Palbo_Soto")]
    for g, col, ls, lab in (("scrb", S.INK, "-", "control"), ("siCDKN1B", S.RUST, (0, (5, 2)), "p27 depleted")):
        d = T[T.genetic_condition == g].sort_values("day")
        ax.plot(d.day, d.published_mean_approx / float(d[d.day == 3].published_mean_approx.iloc[0]), color=col, ls=ls, lw=1.3, marker="o", ms=3, label=lab)
    ax.axvline(3, color=S.GREY, lw=0.5, ls=(0, (1, 2)))
    ax.set_ylim(0.5, 2.05)
    ax.text(3.15, 0.53, "withdrawal", fontsize=6.8, color=S.GREY, va="bottom")
    ax.set_xlabel("Day")
    ax.set_ylabel("Growth relative to day 3")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, fontsize=7.5,
              handlelength=1.6, columnspacing=0.7, borderaxespad=0)


def panel_ratio(ax):
    K = pd.read_csv(PANC_TAB + "/combination_key_results.csv")
    labels = {"F": "low dose", "G": "middle dose", "H": "high dose", "I": "Murine-cell\nconfluence", "J": "Murine-cell\nviability"}
    for i, r in K.iterrows():
        v = r.depleted_vs_intact_recovery_fold_ratio
        # This is a descriptive ratio, not a test of agreement with the mechanism.
        col = S.INK
        ax.plot([1, v], [i, i], color=col, lw=0.9)
        ax.plot(v, i, marker="o", ms=5, mfc=col if r.panel in "FGH" else "white", mec=col, mew=1.0, zorder=3)
        ax.text(v + 0.07 if v > 1 else 1.07, i, f"{v:.2f}", ha="left", va="center", fontsize=7)
    ax.axvline(1, color=S.GREY, lw=0.5, ls=(0, (3, 2)))
    ax.set_yticks(range(len(K)))
    ax.set_yticklabels([labels[p] for p in K.panel], fontsize=7.5)
    ax.set_ylim(len(K) - 0.5, -0.5)
    ax.set_xlim(0.5, 1.8)
    ax.set_xlabel("Fold-change ratio,\np27 loss / control")


def panel_cells(ax, gse, phen, groups, title):
    M = pd.read_csv(TAB + "/M1_cell_line_margins.csv").query("dataset == @gse")
    xs = []
    for j, p in enumerate(phen):
        for k, (g, col) in enumerate(zip(groups, (S.GREY, S.RUST))):
            v = M[(M.phenotype == p) & (M.group == g)].margin.values
            x = j * 2.6 + k
            xs.append((x, g))
            jitter = np.linspace(-0.12, 0.12, len(v))
            ax.scatter(x + jitter, v, s=16, facecolor="white" if k == 0 else col, edgecolor=col, lw=1.0, zorder=3)
            ax.plot([x - 0.25, x + 0.25], [v.mean()] * 2, color=S.INK, lw=1.1)
        a = M[(M.phenotype == p) & (M.group == groups[0])].margin.mean()
        b = M[(M.phenotype == p) & (M.group == groups[1])].margin.mean()
        ax.text(j * 2.6 + 0.5, 4.2, f"{p}\nshift {b - a:+.2f}".replace("-", "\u2212"), ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color=S.GREY, lw=0.5, ls=(0, (3, 2)))
    ax.set_xticks([x for x, _ in xs])
    ax.set_xticklabels([g for _, g in xs], fontsize=7, rotation=30, ha="right")
    ax.set_xlim(-0.6, max(x for x, _ in xs) + 0.6)
    ax.set_ylim(-6, 6.6)
    ax.set_ylabel("Margin (ON minus OFF)")


def panel_cptac(ax):
    C = pd.read_csv(TAB + "/M2_cptac_contexts.csv")
    R = pd.read_csv(TAB + "/M3_cptac_readouts.csv").set_index("readout").loc["prot_CDKN1B"]
    rng = np.random.default_rng(2)
    for k, (flag, lab, col) in enumerate(((False, "other\ncontexts", S.GREY), (True, "proliferative\ncontexts", S.RUST))):
        v = C.loc[C.proliferative_context == flag, "prot_CDKN1B"].dropna().values
        ax.scatter(k + rng.uniform(-0.12, 0.12, len(v)), v, s=7, color=col, alpha=0.55, lw=0)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([k + 0.25, k + 0.25], [q1, q3], color=col, lw=2.4, solid_capstyle="butt")
        ax.plot([k - 0.2, k + 0.32], [med, med], color=S.INK, lw=1.1)
    ax.set_xticks([0, 1])
    n0 = int(C.loc[~C.proliferative_context, "prot_CDKN1B"].notna().sum()); n1 = int(C.loc[C.proliferative_context, "prot_CDKN1B"].notna().sum())
    ax.set_xticklabels([f"other\ncontexts (n = {n0})", f"proliferative\ncontexts (n = {n1})"], fontsize=7.3)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-2.6, 2.1)
    ax.set_ylabel("p27 protein (CPTAC)")


if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 7.1))
    w, h = 0.245, 0.215
    xs = [0.075, 0.405, 0.735]
    r1, r2, r3 = 0.735, 0.395, 0.065
    ax = fig.add_axes([xs[0], r1, w, h]); panel_cells(ax, "GSE236228", ["Proliferation", "Apoptosis"], ["control", "penfluridol"], "Ishikawa cells"); S.panel(ax, "a", dx=-0.2, dy=1.13)
    ax = fig.add_axes([xs[1] + 0.06, r1, w - 0.12, h]); panel_cells(ax, "GSE61481", ["Proliferation"], ["aFGF", "aFGF + AZD4547"], "AN3-CA cells"); S.panel(ax, "b", dx=-0.4, dy=1.13)
    ax = fig.add_axes([xs[2], r1, w, h]); panel_cptac(ax); S.panel(ax, "c", dx=-0.2, dy=1.13)
    ax = fig.add_axes([xs[0], r2, w, h]); panel_model(ax); S.panel(ax, "d", dx=-0.2, dy=1.05)
    ax = fig.add_axes([xs[1], r2, w, h]); panel_duration(ax, "phospho_Rb", "Phospho-Rb positive (%)", 85); S.panel(ax, "e", dx=-0.2, dy=1.05)
    ax.legend(loc="lower left", fontsize=7, handlelength=1.6, labelspacing=0.2, bbox_to_anchor=(-0.02, -0.02))
    ax = fig.add_axes([xs[2], r2, w, h]); panel_duration(ax, "EdU", "EdU positive (%)", 28); S.panel(ax, "f", dx=-0.2, dy=1.05)
    ax = fig.add_axes([xs[0] + 0.045, r3, w - 0.045, h]); panel_depletion(ax); S.panel(ax, "g", dx=-0.52, dy=1.05)
    ax = fig.add_axes([xs[1], r3, w, h]); panel_traj(ax); S.panel(ax, "h", dx=-0.2, dy=1.05)
    ax = fig.add_axes([xs[2] + 0.07, r3, w - 0.065, h]); panel_ratio(ax); S.panel(ax, "i", dx=-0.62, dy=1.05)
    S.save(fig, FIGURES, "Figure7")
