"""Figure 2. The model and its input contexts.
(a) the full network: 37 nodes and 168 signed interactions drawn with horizontal and vertical connectors; (b) the 16 input contexts with the number of non-input nodes fixed by
the inputs and the number of attractors; (c) the seven free nodes of C0101 (orthogonal edges) and the states of its
two fixed points."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

import ecbn
import netdraw as N
import style as S
from paths import FIGURES, TAB

S.apply()
def panel_contexts(ax_code, ax_fix, ax_h):
    D = pd.read_csv(TAB + "/K1_context_landscape.csv")
    D = D.sort_values(["n_proliferative_attractors", "nodes_fixed_by_inputs"], ascending=[False, False]).reset_index(drop=True)
    D["prolif"] = D.n_proliferative_attractors > 0
    D = pd.concat([D[D.prolif].sort_values("nodes_fixed_by_inputs", ascending=False), D[~D.prolif].sort_values("nodes_fixed_by_inputs", ascending=False)]).reset_index(drop=True)
    y = np.arange(len(D))
    for i, ctx in enumerate(D.context):
        for j, b in enumerate(ctx[1:]):
            ax_code.add_patch(Rectangle((j - 0.4, i - 0.36), 0.8, 0.72, fc=S.INK if b == "1" else "white", ec=S.INK, lw=0.4))
    ax_code.set_xlim(-0.6, 3.6)
    ax_code.set_xticks(range(4))
    ax_code.set_xticklabels(ecbn.INPUTS, rotation=90, fontsize=6.5)
    ax_code.xaxis.tick_top()
    ax_code.tick_params(length=0, pad=1.5)
    ax_code.set_yticks(y)
    ax_code.set_yticklabels(D.context, fontsize=6.5)
    for s in ax_code.spines.values():
        s.set_visible(False)
    col = [S.BLUE if p else S.LIGHT for p in D.prolif]
    ax_fix.barh(y, D.nodes_fixed_by_inputs, color=col, height=0.72)
    ax_fix.set_xlim(0, 33)
    ax_fix.set_xticks([0, 11, 22, 33])
    ax_fix.set_xlabel("Nodes fixed by\nthe inputs (of 33)", fontsize=8)
    ax_h.barh(y, D.n_attractors, color=col, height=0.72)
    ax_h.set_xlim(0, 4.4)
    ax_h.set_xticks([0, 2, 4])
    ax_h.set_xlabel("Attractors", fontsize=8)
    for a in (ax_fix, ax_h):
        a.set_yticks([])
        a.spines["left"].set_visible(False)
    for a in (ax_code, ax_fix, ax_h):
        a.set_ylim(len(D) - 0.5, -0.5)
    k = int(D.prolif.sum())
    for a in (ax_fix, ax_h):
        a.axhline(k - 0.5, color=S.GREY, lw=0.5, ls=(0, (3, 2)))
    from matplotlib.patches import Patch

    ax_fix.legend(handles=[Patch(color=S.BLUE, label="context with a proliferative attractor"), Patch(color=S.LIGHT, label="context without")],
                  loc="upper left", bbox_to_anchor=(-0.85, -0.29), ncol=1, fontsize=7, handlelength=1.0, handletextpad=0.4, labelspacing=0.2)


GROW = {"AKT2": 1, "CCND2": 1, "CCND3": 1, "CDKN1B": 0, "IRF1": 0, "MYC": 1, "Proliferation": 1}
ARREST = {k: 1 - v for k, v in GROW.items()}

def state_table(ax):
    A3 = pd.read_csv(TAB + "/A3_context_attractor_probabilities.csv").set_index("attractor").MaBoSS_probability_t1000
    nodes = ["MYC", "AKT2", "CCND2", "CCND3", "CDKN1B", "IRF1", "Proliferation"]
    rows = [("proliferative", GROW, A3["C0101_A1"]), ("arrested", ARREST, A3["C0101_A2"])]
    for r, (name, st, w) in enumerate(rows):
        for c, n in enumerate(nodes):
            ax.add_patch(Rectangle((c - 0.4, r - 0.36), 0.8, 0.72, fc=S.INK if st[n] else "white", ec=S.INK, lw=0.5))
        ax.text(-0.75, r, name, ha="right", va="center", fontsize=8)
        ax.text(len(nodes) + 0.75, r, f"{w:.2f}", ha="center", va="center", fontsize=8)
    for c, n in enumerate(nodes):
        ax.text(c, -0.6, S.label(n), ha="left", va="bottom", rotation=50, fontsize=7.5)
    ax.text(len(nodes) + 0.75, -0.6, "basin\nweight", ha="center", va="bottom", fontsize=7.5)
    ax.text(-0.75, -0.75, "fixed point", ha="right", va="bottom", fontsize=7.5, style="italic")
    ax.set_xlim(-3.6, len(nodes) + 2.0)
    ax.set_ylim(1.7, -2.3)
    ax.set_aspect("equal")
    ax.axis("off")


SWITCH = {"IRF1", "MYC", "AKT2", "CDKN1B", "CCND2", "CCND3", "Proliferation"}


def panel_network(fig, y0, W=7.2):
    """All 168 interactions; the seven switch nodes of context C0101 (panel c) are shaded."""
    from matplotlib.lines import Line2D

    boxes, bands, H = N.full_layout(W, py=0.27, head=0.30, foot=0.04)
    ax = N.net_axes(fig, 0, y0, W, H)
    N.draw_bands(ax, bands, fs=6.8)
    N.draw_edges(ax, boxes, N.signed_edges(), lw=0.40, head=0.024, lane_sep=0.014, track_sep=0.014)

    def style(n):
        fc, ec, tc, lw, ls = N.default_node_style(n)
        return (S.LIGHT if n in SWITCH else fc), ec, tc, lw, ls

    N.draw_nodes(ax, boxes, style)
    ax.legend(handles=[Line2D([], [], color=N.ACT, lw=0.9, marker=">", ms=3.5, markevery=[1], label="activation"),
                       Line2D([], [], color=N.INH, lw=0.9, marker="|", ms=5, markevery=[1], label="inhibition")],
              loc="upper center", bbox_to_anchor=(0.5, -0.005), ncol=2, fontsize=7.5, handlelength=2.2,
              columnspacing=1.6, borderaxespad=0.0)
    return ax, H


if __name__ == "__main__":
    W, FH = 7.2, 7.40
    fig = plt.figure(figsize=(W, FH))
    ax_m, H = panel_network(fig, FH - 0.02 - 3.17)
    S.panel(ax_m, "a", dx=0.03, dy=1.0)
    inch = lambda x, y, w, h: fig.add_axes([x / W, y / FH, w / W, h / FH])
    ax_c = inch(0.62, 0.85, 0.42, 2.30)
    ax_f = inch(1.14, 0.85, 0.78, 2.30)
    ax_h = inch(2.07, 0.85, 0.78, 2.30)
    panel_contexts(ax_c, ax_f, ax_h)
    S.panel(ax_c, "b", dx=-1.25, dy=1.20)
    S.attach("b", ax_f, ax_h)
    ax1 = N.net_axes(fig, 3.55, 1.45, 3.55, 1.85)
    N.draw_switch(ax1, 3.55, 1.85)
    S.panel(ax1, "c", dx=0.0, dy=1.0)
    ax2 = inch(3.75, 0.05, 3.35, 1.20)
    state_table(ax2)
    S.attach("c", ax2)
    S.save(fig, FIGURES, "Figure2")
