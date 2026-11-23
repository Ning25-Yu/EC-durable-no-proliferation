"""Figure 4. Suppression during exposure, commitment to arrest, schedules and rate heterogeneity.
(a) Proliferation OFF and commitment to arrest during an L = 8 exposure and after withdrawal; (b) concurrent and
sequential suppression of MYC and cyclin D2 at a matched budget; (c) durable arrest under 200 random rate sets."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style as S
from paths import FIGURES, TAB

S.apply()
T = pd.read_csv(TAB + "/K4_commitment_timecourse.csv")
CASES = [("C0101_A1 | CDKN1B=1", "p27 ON"), ("C0101_A1 | CDKN1A=1", "p21 ON"), ("C0101_A1 | MYC=0", "MYC OFF")]


def panel_time(ax, case, title, first):
    d = T[T.case == case].sort_values("time")
    ax.axvspan(0, 8, color=S.PALE, lw=0)
    ax.plot(d.time, d.proliferation_off, color=S.RUST, lw=1.3, label="Proliferation OFF")
    ax.plot(d.time, d.committed_to_arrest, color=S.BLUE, lw=1.3, ls=(0, (5, 2)), label="committed to arrest")
    ax.set_xlim(0, 30)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xticks([0, 8, 16, 24])
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.set_xlabel("Time")
    ax.text(4, 0.03, "exposure", ha="center", va="bottom", fontsize=7, color=S.GREY)
    if first:
        ax.set_ylabel("Probability")
    else:
        ax.set_yticklabels([])


def panel_schedules(ax):
    D = pd.read_csv(TAB + "/K6_schedules.csv")
    star = pd.read_csv(TAB + "/D4c_time_window_crossover_Lstar.csv").iloc[0]
    ax.plot(D.L, D.concurrent, color=S.BLUE, lw=1.3, label="concurrent")
    ax.plot(D.L, D.MYC_first, color=S.INK, lw=1.1, ls=(0, (5, 2)), label="MYC first")
    ax.plot(D.L, D.CCND2_first, color=S.GREY, lw=1.1, ls=(0, (1.5, 1.5)), label="CCND2 first")
    ax.plot(star.L_star, star.arrest_at_L_star, marker="o", ms=4, mfc="white", mec=S.INK, mew=0.8, zorder=5)
    ax.annotate(f"crossover at {star.L_star:.2f}", (star.L_star, star.arrest_at_L_star), xytext=(6, -9), textcoords="offset points", fontsize=7.5)
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
    ax.set_xlim(0.25, 8)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("Length of each phase")
    ax.set_ylabel("Durable arrest probability")
    ax.legend(loc="upper left", fontsize=7.5, handlelength=2.2, labelspacing=0.25)


def panel_rates(ax):
    R = pd.read_csv(TAB + "/K7_rate_draws.csv")
    series = [
        ("p27_pulse_L4", "p27 ON\nexposure 4", S.GREY),
        ("p27_pulse_L8", "p27 ON\nexposure 8", S.BLUE),
        ("concurrent_L8", "concurrent\nexposure 8", S.BLUE),
        ("MYC_then_CCND2_L8", "MYC first\nexposure 8", S.INK),
        ("CCND2_then_MYC_L8", "CCND2 first\nexposure 8", S.GREY),
    ]
    rng = np.random.default_rng(1)
    for i, (c, lab, col) in enumerate(series):
        v = R[c].values
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=2.2, color=col, alpha=0.45, lw=0, zorder=2)
        q = np.quantile(v, [0.25, 0.5, 0.75])
        ax.plot([i - 0.3, i + 0.3], [q[1], q[1]], color=S.INK, lw=1.2, zorder=3)
        ax.plot([i, i], [q[0], q[2]], color=S.INK, lw=0.7, zorder=3)
    ax.set_xticks(range(len(series)))
    ax.set_xticklabels([s[1] for s in series], fontsize=7)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlim(-0.6, len(series) - 0.4)
    ax.set_ylabel("Durable arrest probability")


if __name__ == "__main__":
    fig = plt.figure(figsize=(7.2, 4.9))
    for k, (case, title) in enumerate(CASES):
        ax = fig.add_axes([0.075 + k * 0.315, 0.60, 0.275, 0.31])
        panel_time(ax, case, title, k == 0)
        if k == 0:
            S.panel(ax, "a", dx=-0.17, dy=1.06)
        else:
            S.attach("a", ax)
        if k == 2:
            ax.legend(loc="center right", fontsize=7.5, handlelength=2.2, labelspacing=0.25)
    ax_b = fig.add_axes([0.075, 0.10, 0.36, 0.34])
    panel_schedules(ax_b)
    S.panel(ax_b, "b", dx=-0.13, dy=1.04)
    ax_c = fig.add_axes([0.56, 0.10, 0.43, 0.34])
    panel_rates(ax_c)
    S.panel(ax_c, "c", dx=-0.11, dy=1.04)
    S.save(fig, FIGURES, "Figure4")
