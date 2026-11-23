"""Shared figure style: restrained palette, serif mathematical type, thin open axes."""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

INK = "#1a1a1a"
GREY = "#7d828a"
LIGHT = "#d5d8dc"
PALE = "#f0f1f3"
BLUE = "#2b5c8a"  # activation, low-risk group, single series
RUST = "#b0522d"  # inhibition, high-risk group, second series
SAND = "#c9a66b"  # third series, used sparingly

PROB = LinearSegmentedColormap.from_list("prob", ["#ffffff", "#c9d6e3", "#6f93b5", BLUE, "#16324d"])

NAME = {"CDKN1A": "p21", "CDKN1B": "p27", "CDKN2A": "p16", "Proliferation": "Proliferation"}
MODULES = [
    ("Inputs", ["KLRB1", "IRF5", "STAT6", "TSC22D3"]),
    ("Interferon", ["IFNG", "IRF1", "CXCL10", "LTA", "SOCS1", "SOCS3"]),
    ("Inflammatory", ["NFKB1", "MAPK14", "IL6", "STAT3", "CCL2", "CCL4", "CCL5", "CXCL8", "TNFRSF11B"]),
    ("Growth and cell cycle", ["IGF1", "AKT2", "MYC", "CCND2", "CCND3", "CDKN1A", "CDKN1B", "CDKN2A", "KLF4", "SPDEF"]),
    ("Effectors", ["TFF3", "MMP2"]),
    ("Phenotypes", ["Proliferation", "Apoptosis", "EMT", "Invasion", "Migration", "Metastasis"]),
]


def sci(p):
    """1.4e-05 -> '1.4 × 10^-5', with only the exponent raised."""
    m, e = f"{p:.1e}".split("e")
    return f"{m} \u00d7 10$^{{{int(e)}}}$"


def label(node):
    return NAME.get(node, node)


def pulse_label(p):
    """'CDKN1B=1' -> 'p27 ON'; 'MYC=0 + CCND2=0' -> 'MYC + CCND2 OFF' style labels."""
    parts = [x.strip() for x in p.split("+")]
    names = [label(x.split("=")[0]) for x in parts]
    vals = [x.split("=")[1] for x in parts]
    if len(set(vals)) == 1:
        return " + ".join(names) + (" ON" if vals[0] == "1" else " OFF")
    return " + ".join(f"{n} {'ON' if v == '1' else 'OFF'}" for n, v in zip(names, vals))


def apply():
    mpl.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "lines.linewidth": 1.2,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def panel(ax, letter, dx=-0.02, dy=1.04):
    t = ax.text(dx, dy, f"({letter})", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")
    fig = ax.figure
    fig.__dict__.setdefault("_panels", {}).setdefault(letter, []).append(ax)
    fig.__dict__.setdefault("_panel_labels", []).append(t)


def attach(letter, *axes):
    """Add further axes (a colour bar, a table, the other members of a row) to panel `letter`."""
    fig = axes[0].figure
    fig.__dict__.setdefault("_panels", {}).setdefault(letter, []).extend(axes)


def _ink_bbox(fig, bb, dpi=300):
    """Shrink a bounding box (inches) to the non-white pixels drawn inside it."""
    import io

    import numpy as np
    from matplotlib.transforms import Bbox
    from PIL import Image

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches=bb, pad_inches=0)
    a = np.asarray(Image.open(buf).convert("L")) < 250
    rows, cols = np.where(a.any(axis=1))[0], np.where(a.any(axis=0))[0]
    if not len(rows):
        return bb
    h = a.shape[0]
    return Bbox.from_extents(bb.x0 + cols[0] / dpi, bb.y0 + (h - 1 - rows[-1]) / dpi,
                             bb.x0 + (cols[-1] + 1) / dpi, bb.y0 + (h - rows[0]) / dpi)


def save_panels(fig, folder, prefix, pad=0.02):
    """Save every labelled panel as its own PDF, cropped to its content, without the panel letter."""
    from matplotlib.transforms import Bbox

    os.makedirs(folder, exist_ok=True)
    groups = fig.__dict__.get("_panels", {})
    loose = [ax for ax in fig.axes if not any(ax in g for g in groups.values())]
    assert not loose, f"{prefix}: axes not assigned to a panel"
    for t in fig.__dict__.get("_panel_labels", []):
        t.set_visible(False)
    shown = {ax: ax.get_visible() for ax in fig.axes}
    for letter, group in sorted(groups.items()):
        for ax in fig.axes:
            ax.set_visible(ax in group)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        bb = Bbox.union([ax.get_tightbbox(r) for ax in group]).transformed(fig.dpi_scale_trans.inverted())
        bb = _ink_bbox(fig, bb)  # drawing axes (axis off) are larger than what is drawn in them
        bb = Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        fig.savefig(os.path.join(folder, f"{prefix}{letter}.pdf"), bbox_inches=bb)
    for ax, v in shown.items():
        ax.set_visible(v)
    for t in fig.__dict__.get("_panel_labels", []):
        t.set_visible(True)


def save(fig, folder, name):
    """Main-text figures (Figure1-Figure7) are saved only as one PDF per panel (fig1a.pdf, fig1b.pdf, ...), which
    the manuscript assembles; supplementary figures are saved whole as PDF and 600 dpi PNG."""
    os.makedirs(folder, exist_ok=True)
    if name[len("Figure"):].isdigit():
        save_panels(fig, folder, "fig" + name[len("Figure"):])
    else:
        fig.savefig(os.path.join(folder, name + ".pdf"), bbox_inches="tight", pad_inches=0.03)
        fig.savefig(os.path.join(folder, name + ".png"), bbox_inches="tight", pad_inches=0.03, dpi=600)
    plt.close(fig)
