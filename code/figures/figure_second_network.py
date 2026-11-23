"""Figure S4. The same analysis on a second cancer model, the reduced T-LGL leukaemia survival network.
(a) the six nodes with horizontal and vertical connectors (filled, ON in the survival fixed point); (b) probability
a(L) that a pulse applied in the survival fixed point leaves the network in apoptosis, closed forms and matrix
exponential; (c) exposure durations L50 and L90 of the five interventions that work."""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import netdraw as N
import style as S
from orthonet import Box
from paths import FIGURES, TAB

S.apply()

EDGES = [("Ceramide", "S1P", True), ("DISC", "FLIP", True), ("S1P", "Fas", True), ("Fas", "Ceramide", False),
         ("S1P", "Ceramide", True), ("Ceramide", "DISC", False), ("Fas", "DISC", False), ("FLIP", "DISC", True),
         ("DISC", "Apoptosis", False)]
GRID = {"S1P": (0, 0), "Ceramide": (1, 0), "DISC": (2, 0), "Apoptosis": (3, 0), "Fas": (0, 1), "FLIP": (2, 1)}
ON = {"S1P", "FLIP"}
STYLE = {
    "DISC ON": (S.INK, "-"),
    "Ceramide ON": (S.INK, (0, (4, 1.6))),
    "S1P OFF": (S.BLUE, "-"),
    "S1P OFF + FLIP OFF": (S.BLUE, (0, (4, 1.6))),
    "Fas ON + FLIP OFF": (S.RUST, "-"),
    "Fas ON": (S.GREY, "-"),
    "FLIP OFF": (S.GREY, "-"),
}


def main():
    D = pd.read_csv(TAB + "/J2_tlgl_duration_response.csv")
    fig = plt.figure(figsize=(7.2, 2.55))
    W, H = 2.95, 1.55
    ax = N.net_axes(fig, 0.05, 0.50, W, H)
    bw, bh = 0.50, 0.19
    px, py = (W - 0.16 - bw) / 3, (H - bh - 0.2) / 1
    boxes = {n: Box(0.08 + bw / 2 + c * px, H - 0.1 - bh / 2 - r * py, bw, bh) for n, (c, r) in GRID.items()}
    N.draw_edges(ax, boxes, EDGES, lw=0.7, head=0.04, lane_sep=0.05, track_sep=0.045, margin=0.11, port_frac=0.72)
    for n, b in boxes.items():
        on = n in ON
        N.node(ax, b, n, S.INK if on else "white", S.INK, "white" if on else S.INK, 6.4, 0.7)
    S.panel(ax, "a", dx=0.0, dy=1.17)

    axb = fig.add_axes([0.505, 0.175, 0.20, 0.70])
    axc = fig.add_axes([0.83, 0.175, 0.16, 0.70])
    for name, (col, ls) in STYLE.items():
        d = D[D.pulse == name]
        axb.plot(d.L, d.arrest, color=col, ls=ls, lw=1.25)
        m = d[(d.L > 0) & (d.L * 2 % 2 == 0) & (d.L <= 8)]
        if d.arrest.max() > 0:
            axb.plot(m.L, m.arrest_expm, "o", ms=2.6, mfc="white", mec=col, mew=0.7)
    axb.set_xlim(0, 8)
    axb.set_ylim(-0.03, 1.03)
    axb.set_xlabel("Exposure duration (model time)")
    axb.set_ylabel("Durable apoptosis probability")
    axb.set_xticks([0, 2, 4, 6, 8])
    J3 = pd.read_csv(TAB + "/J3_tlgl_summary.csv").dropna(subset=["L90"]).sort_values("L90")
    for i, r in enumerate(J3.itertuples()):
        col, ls = STYLE[r.pulse]
        axc.plot([r.L50, r.L90], [i, i], color=S.GREY, lw=0.8, zorder=1)
        axc.plot(r.L50, i, marker="o", ms=5, mfc="white", mec=col, mew=1.0, zorder=3)
        axc.plot(r.L90, i, marker="o", ms=5, mfc=col, mec=col, zorder=3)
    axc.set_yticks(range(len(J3)))
    axc.set_yticklabels([p.replace(" + ", "\n+ ") for p in J3.pulse], fontsize=7)
    axc.set_ylim(len(J3) - 0.5, -0.5)
    axc.set_xlim(-0.2, 4)
    axc.set_xlabel("Exposure duration")
    axc.plot([], [], ls="none", marker="o", ms=5, mfc="white", mec=S.INK, mew=1.0, label="50% apoptosis")
    axc.plot([], [], ls="none", marker="o", ms=5, mfc=S.INK, mec=S.INK, label="90% apoptosis")
    axc.legend(loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=2, fontsize=7, handletextpad=0.1, columnspacing=1.0, frameon=False)
    S.panel(axb, "b")
    S.panel(axc, "c", dx=-0.55)
    handles = [Line2D([], [], color=c, ls=ls, lw=1.25, label=n) for n, (c, ls) in STYLE.items() if n not in ("Fas ON", "FLIP OFF")]
    handles.append(Line2D([], [], color=S.GREY, lw=1.25, label="Fas ON, FLIP OFF (alone)"))
    axb.legend(handles=handles, loc="lower right", fontsize=5.6, handlelength=2.2, frameon=False, labelspacing=0.3,
               bbox_to_anchor=(1.03, 0.03), borderaxespad=0.2)
    edge = [Line2D([], [], color=N.ACT, lw=0.9, label="activation"), Line2D([], [], color=N.INH, lw=0.9, label="inhibition")]
    fig.legend(handles=edge, loc="lower left", fontsize=6.6, handlelength=2.0, frameon=False, ncol=2, columnspacing=1.4,
               bbox_to_anchor=(0.10, 0.06))
    S.save(fig, FIGURES, "FigureS4")


if __name__ == "__main__":
    main()
