"""Figure 1 panels redrawn for the panel-by-panel layout.
fig1a: the authors' study workflow (wording and colours as in their drawing), set as a tall left column so that
       its text matches the other panels in size;
fig1f: high-risk fraction by context, wide enough to sit under panels d and e.
Both are written to results/figures (replacing the versions from figure1_transcriptomic.py) at the size they are shown in the manuscript (fig1a at scale 1, fig1f at 0.94).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mp
import matplotlib.pyplot as plt

import style as S
from figure1_transcriptomic import panel_contexts
from paths import FIGURES

PANELS = FIGURES
LINE, EDGE, TEXT = "#3e535b", "#87979c", "#13252b"
PT = 1 / 72


# ----------------------------------------------------------------------------- fig1a
ROWS = [  # each row: list of boxes (bold title, regular text); widths are fractions of the panel width
    [("TCGA tumour and\nnormal endometrium", "Differential expression\nwith limma\n7,977 differentially\nexpressed genes", 0.54),
     ("ImmPort", "Immune-related\ngene list", 0.40)],
    [("759 immune-related\ndifferentially expressed genes", None, 0.86)],
    [("Univariable Cox screening\nin all 519 tumours", "108 candidate genes", 0.86)],
    [("7:3 split stratified\nby vital status", None, 0.86)],
    [("Training set", "363 tumours", 0.44), ("Test set", "156 tumours for\ninternal assessment", 0.44)],
    [("LASSO Cox regression", "24 genes retained", 0.56)],
    [("Stepwise Cox regression\nKLRB1, IRF5, STAT6\nand TSC22D3", None, 0.58)],
    [("Risk groups from the\nfour-gene score", "Training set median\napplied to both sets", 0.49),
     ("16 Boolean input\ncontexts", "Each gene uses its\ntraining median", 0.47)],
    [("Survival\ncomparison", "Training and\ntest sets", 0.315), ("Hallmark score\ncomparison", "502 tumours", 0.315),
     ("Risk group by\ninput context", "Fisher's exact\ntest", 0.315)],
]
# x centres (fractions); children of each row (row index -> list of (parent box, child boxes))
XC = [[0.27, 0.80], [0.5], [0.5], [0.5], [0.30, 0.77], [0.30], [0.30], [0.245, 0.765], [0.1575, 0.5, 0.8425]]
LINKS = {0: [([0, 1], [0])], 1: [([0], [0])], 2: [([0], [0])], 3: [([0], [0, 1])], 4: [([0], [0])],
         5: [([0], [0])], 6: [([0], [0, 1])], 7: [([0], [0, 1, 2]), ([1], [2])]}
FS_T, FS_R, LEAD = 7.0, 6.6, 1.22


def box_height(title, body):
    n_t = title.count("\n") + 1
    n_b = body.count("\n") + 1 if body else 0
    return (n_t * FS_T + n_b * FS_R) * LEAD + 2 * 3.2 + (1.5 if body else 0)


def fig1a(W=170.0, H=440.0):
    fig = plt.figure(figsize=(W * PT, H * PT))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    hts = [max(box_height(t, b) for t, b, _ in row) for row in ROWS]
    gap = (H - sum(hts) - 2) / (len(ROWS) - 1)
    print("fig1a row gap (pt):", round(gap, 1))
    boxes, y = [], H - 1
    for row, xc, h in zip(ROWS, XC, hts):
        rb = []
        for (title, body, wf), x in zip(row, xc):
            w = wf * W
            bold_box = title.startswith("Stepwise")
            ax.add_patch(mp.FancyBboxPatch((x * W - w / 2, y - h), w, h, boxstyle="round,pad=0,rounding_size=3",
                                           fc="white", ec=LINE if bold_box else EDGE, lw=1.0 if bold_box else 0.6))
            n_t = title.count("\n") + 1
            n_b = body.count("\n") + 1 if body else 0
            block = (n_t * FS_T + n_b * FS_R) * LEAD + (1.5 if body else 0)
            yt = y - h / 2 + block / 2
            ax.text(x * W, yt, title, ha="center", va="top", fontsize=FS_T, fontweight="bold", color=TEXT,
                    linespacing=LEAD, family="Liberation Sans")
            if body:
                ax.text(x * W, yt - n_t * FS_T * LEAD - 1.5, body, ha="center", va="top", fontsize=FS_R, color=TEXT,
                        linespacing=LEAD, family="Liberation Sans")
            rb.append(dict(x=x * W, top=y, bot=y - h))
        boxes.append(rb)
        y -= h + gap
    head = dict(arrowstyle="-|>,head_length=0.32,head_width=0.16", color=LINE, lw=0.6, shrinkA=0, shrinkB=0.3)
    for r, links in LINKS.items():
        for parents, children in links:
            P = [boxes[r][i] for i in parents]
            C = [boxes[r + 1][i] for i in children]
            ybus = (min(p["bot"] for p in P) + max(c["top"] for c in C)) / 2
            xs = [p["x"] for p in P] + [c["x"] for c in C]
            for p in P:
                ax.plot([p["x"], p["x"]], [p["bot"], ybus], color=LINE, lw=0.6, solid_capstyle="butt")
            if len(set(round(v, 2) for v in xs)) > 1:
                ax.plot([min(xs), max(xs)], [ybus, ybus], color=LINE, lw=0.6, solid_capstyle="butt")
            for c in C:
                ax.annotate("", xy=(c["x"], c["top"]), xytext=(c["x"], ybus), arrowprops=head)
    os.makedirs(PANELS, exist_ok=True)
    fig.savefig(os.path.join(PANELS, "fig1a.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------- fig1f
def fig1f():
    S.apply()
    fig = plt.figure(figsize=(4.55, 1.62))
    ax = fig.add_axes([0.21, 0.30, 0.76, 0.56])
    panel_contexts(ax)
    S.panel(ax, "f")
    S.save_panels(fig, PANELS, "fig1")
    plt.close(fig)


if __name__ == "__main__":
    fig1a()
    fig1f()
