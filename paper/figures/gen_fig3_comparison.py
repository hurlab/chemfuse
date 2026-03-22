"""
Generate Fig. 3 – ChemFuse feature comparison matrix.

Publication-quality heatmap, saved as both SVG and PDF.
Run with: .venv/bin/python paper/figures/gen_fig3_comparison.py
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

matplotlib.use("Agg")  # headless / no display needed

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
tools = [
    "ChemFuse",
    "PubChemPy",
    "RDKit",
    "BioServices",
    "SwissADME",
    "ADMETlab 3.0",
    "admet-ai",
    r"Schr$\"{o}$dinger Suite",  # rendered inline; raw label used for file
]

features = [
    "Multi-DB\nSearch",
    "Descriptors",
    "ADMET",
    "Drug-\nlikeness",
    "Fingerprints",
    "Clustering",
    "Python\nAPI",
    "R\nInterface",
    "CLI",
    "Web\nUI",
    "Open\nSource",
    "Free",
]

data = np.array(
    [
        [1,    1,    1,    1,    1,    1,    1,    0.5,  1,    1,    1,    1   ],  # ChemFuse
        [0.5,  0,    0,    0,    0,    0,    1,    0,    0,    0,    1,    1   ],  # PubChemPy
        [0,    1,    0,    1,    1,    1,    1,    0,    0,    0,    1,    1   ],  # RDKit
        [1,    0,    0,    0,    0,    0,    1,    0,    0,    0,    1,    1   ],  # BioServices
        [0,    0.5,  1,    1,    0,    0,    0,    0,    0,    1,    0,    1   ],  # SwissADME
        [0,    0,    1,    0.5,  0,    0,    0.5,  0,    0,    1,    0,    1   ],  # ADMETlab 3.0
        [0,    0,    1,    0,    0,    0,    1,    0,    1,    1,    1,    1   ],  # admet-ai
        [0.5,  1,    1,    1,    1,    1,    1,    0.5,  1,    1,    0,    0   ],  # Schrodinger
    ]
)

n_tools, n_feat = data.shape

# ---------------------------------------------------------------------------
# Color map: 0 -> white/light-gray | 0.5 -> light green | 1 -> dark green
# ---------------------------------------------------------------------------
colors_list = ["#F0F0F0", "#A8D5A2", "#2D7D32"]   # none / partial / full
cmap = ListedColormap(colors_list)
bounds = [-0.25, 0.25, 0.75, 1.25]
norm = BoundaryNorm(bounds, cmap.N)

# ---------------------------------------------------------------------------
# Figure layout
# 180 mm x 120 mm in inches (1 inch = 25.4 mm)
# ---------------------------------------------------------------------------
fig_w = 180 / 25.4   # ~7.09 in
fig_h = 120 / 25.4   # ~4.72 in

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
fig.patch.set_facecolor("white")

im = ax.imshow(data, cmap=cmap, norm=norm, aspect="auto")

# ---------------------------------------------------------------------------
# Cell annotations
# ---------------------------------------------------------------------------
annotation_font = dict(
    fontsize=8,
    fontfamily="DejaVu Sans",
    va="center",
    ha="center",
)

for row in range(n_tools):
    for col in range(n_feat):
        val = data[row, col]
        if val == 1.0:
            symbol = "\u2713"           # checkmark U+2713
            color = "white"
            weight = "bold"
        elif val == 0.5:
            symbol = "P"
            color = "#1B5E20"           # dark text on light green
            weight = "bold"
        else:
            symbol = "\u2013"           # en-dash
            color = "#9E9E9E"           # gray text on light-gray cell
            weight = "normal"

        ax.text(
            col, row, symbol,
            color=color,
            fontweight=weight,
            **annotation_font,
        )

# ---------------------------------------------------------------------------
# Axes – ticks and labels
# ---------------------------------------------------------------------------
ax.set_xticks(np.arange(n_feat))
ax.set_yticks(np.arange(n_tools))

# Column labels – rotated 45 degrees, right-aligned
ax.set_xticklabels(
    features,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=7.5,
)

# Row labels – plain strings (no LaTeX needed for Schrodinger here)
tool_labels_plain = [
    "ChemFuse",
    "PubChemPy",
    "RDKit",
    "BioServices",
    "SwissADME",
    "ADMETlab 3.0",
    "admet-ai",
    "Schr\u00f6dinger Suite",   # unicode ö
]
ax.set_yticklabels(tool_labels_plain, fontsize=8)

# Bold the ChemFuse row label (index 0)
ytick_labels = ax.get_yticklabels()
ytick_labels[0].set_fontweight("bold")
ytick_labels[0].set_fontsize(8.5)

# Move x-axis ticks/labels to top for a cleaner matrix look
ax.xaxis.set_label_position("top")
ax.xaxis.tick_top()
# Rotate top labels (must redo after tick_top)
ax.set_xticklabels(
    features,
    rotation=45,
    ha="left",
    rotation_mode="anchor",
    fontsize=7.5,
)

# ---------------------------------------------------------------------------
# Grid lines – thin white lines between cells
# ---------------------------------------------------------------------------
ax.set_xticks(np.arange(n_feat) - 0.5, minor=True)
ax.set_yticks(np.arange(n_tools) - 0.5, minor=True)
ax.grid(which="minor", color="white", linewidth=1.2)
ax.tick_params(which="minor", bottom=False, left=False, top=False)

# Remove outer frame spines for clean academic look
for spine in ax.spines.values():
    spine.set_visible(False)

# Remove major tick marks
ax.tick_params(which="major", length=0)

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
patch_full    = mpatches.Patch(facecolor="#2D7D32", edgecolor="#CCCCCC",
                                linewidth=0.5, label="Full support  (\u2713)")
patch_partial = mpatches.Patch(facecolor="#A8D5A2", edgecolor="#CCCCCC",
                                linewidth=0.5, label="Partial  (P)")
patch_none    = mpatches.Patch(facecolor="#F0F0F0", edgecolor="#CCCCCC",
                                linewidth=0.5, label="Not available  (\u2013)")

legend = ax.legend(
    handles=[patch_full, patch_partial, patch_none],
    loc="lower right",
    bbox_to_anchor=(1.0, -0.32),
    ncol=3,
    fontsize=7.5,
    frameon=True,
    framealpha=0.95,
    edgecolor="#CCCCCC",
    title="Support level",
    title_fontsize=7.5,
)

# ---------------------------------------------------------------------------
# Tight layout and save
# ---------------------------------------------------------------------------
fig.tight_layout(rect=[0, 0.06, 1, 0.97])

out_dir = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures"
svg_path = f"{out_dir}/fig3_comparison.svg"
pdf_path = f"{out_dir}/fig3_comparison.pdf"

fig.savefig(svg_path, format="svg", dpi=300, bbox_inches="tight",
            facecolor="white")
fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight",
            facecolor="white")

print(f"Saved SVG: {svg_path}")
print(f"Saved PDF: {pdf_path}")
plt.close(fig)
