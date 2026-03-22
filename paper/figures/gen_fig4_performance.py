"""
Generate Figure 4: Performance comparison chart for ChemFuse paper.
Dual-panel figure:
  Panel A - Batch Screening Time (ChemFuse vs Manual workflow)
  Panel B - Cache Impact (Cold vs Warm cache)
"""

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless generation

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

compounds = [100, 500, 1000]
x = np.array(compounds)

# Panel A
chemfuse_times = np.array([0.75, 3.0, 5.7])
manual_times   = np.array([3.3,  16.0, 35.0])
speedup_factors = ["4.4x", "5.3x", "6.2x"]

# Panel B
cold_cache = np.array([0.75, 3.0, 5.7])
warm_cache = np.array([0.12, 0.5, 0.75])
cache_speedup = [f"{c/w:.1f}x" for c, w in zip(cold_cache, warm_cache)]

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------

FONT_FAMILY   = "DejaVu Sans"
LABEL_SIZE    = 9
TICK_SIZE     = 8
ANNOT_SIZE    = 8
PANEL_LABEL   = 10

COLOR_CHEMFUSE  = "#2166AC"   # blue
COLOR_MANUAL    = "#D6604D"   # red/orange
COLOR_COLD      = "#2166AC"   # darker blue
COLOR_WARM      = "#92C5DE"   # lighter blue
COLOR_GRID      = "#D0D0D0"

plt.rcParams.update({
    "font.family":        FONT_FAMILY,
    "font.size":          LABEL_SIZE,
    "axes.labelsize":     LABEL_SIZE,
    "axes.titlesize":     LABEL_SIZE,
    "xtick.labelsize":    TICK_SIZE,
    "ytick.labelsize":    TICK_SIZE,
    "legend.fontsize":    TICK_SIZE,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "lines.linewidth":    1.6,
    "lines.markersize":   6,
})

# Figure size: ~180mm x 80mm converted to inches (1 inch = 25.4 mm)
fig_width  = 180 / 25.4   # ~7.09 in
fig_height =  80 / 25.4   # ~3.15 in

fig, (ax_a, ax_b) = plt.subplots(
    1, 2,
    figsize=(fig_width, fig_height),
    constrained_layout=True,
)

# ---------------------------------------------------------------------------
# Panel A: Batch Screening Time
# ---------------------------------------------------------------------------

# Shaded region between the two lines
ax_a.fill_between(
    x, chemfuse_times, manual_times,
    alpha=0.12, color=COLOR_CHEMFUSE, label=None
)

# Manual workflow line
ax_a.plot(
    x, manual_times,
    color=COLOR_MANUAL, marker="s", linestyle="-",
    label="Manual workflow", zorder=3,
)

# ChemFuse line
ax_a.plot(
    x, chemfuse_times,
    color=COLOR_CHEMFUSE, marker="o", linestyle="-",
    label="ChemFuse", zorder=3,
)

# Speedup annotations at midpoint of gap
for xi, ct, mt, sf in zip(x, chemfuse_times, manual_times, speedup_factors):
    mid_y = (ct + mt) / 2.0
    ax_a.annotate(
        sf,
        xy=(xi, mid_y),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=ANNOT_SIZE,
        color="#444444",
        va="center",
        ha="left",
    )

ax_a.set_xlabel("Number of compounds")
ax_a.set_ylabel("Time (minutes)")
ax_a.set_xticks(x)
ax_a.set_xticklabels([str(c) for c in compounds])
ax_a.set_ylim(bottom=0)
ax_a.yaxis.grid(True, color=COLOR_GRID, linestyle="--", linewidth=0.6, zorder=0)
ax_a.set_axisbelow(True)
ax_a.legend(loc="upper left", framealpha=0.85, edgecolor="#CCCCCC")

# Disclaimer note
ax_a.text(
    0.98, 0.03,
    "Approx. values; actual times\ndepend on network conditions",
    transform=ax_a.transAxes,
    fontsize=6.5,
    color="#888888",
    ha="right", va="bottom",
    linespacing=1.4,
)

# Panel label (A)
ax_a.text(
    -0.14, 1.04, "(A)",
    transform=ax_a.transAxes,
    fontsize=PANEL_LABEL, fontweight="bold",
    va="top", ha="left",
)

# ---------------------------------------------------------------------------
# Panel B: Cache Impact
# ---------------------------------------------------------------------------

n_groups  = len(compounds)
bar_width = 0.32
offsets   = np.array([-bar_width / 2, bar_width / 2])

bar_positions_cold = np.arange(n_groups) + offsets[0]
bar_positions_warm = np.arange(n_groups) + offsets[1]

bars_cold = ax_b.bar(
    bar_positions_cold, cold_cache,
    width=bar_width, color=COLOR_COLD,
    label="Cold cache", zorder=3,
)
bars_warm = ax_b.bar(
    bar_positions_warm, warm_cache,
    width=bar_width, color=COLOR_WARM,
    label="Warm cache", zorder=3,
    edgecolor=COLOR_COLD, linewidth=0.5,
)

# Speedup annotations above warm cache bars
for pos, wt, sf in zip(bar_positions_warm, warm_cache, cache_speedup):
    ax_b.text(
        pos, wt + 0.25,
        sf,
        ha="center", va="bottom",
        fontsize=ANNOT_SIZE,
        color="#444444",
    )

ax_b.set_xlabel("Number of compounds")
ax_b.set_ylabel("Time (minutes)")
ax_b.set_xticks(np.arange(n_groups))
ax_b.set_xticklabels([str(c) for c in compounds])
ax_b.set_ylim(bottom=0)
ax_b.yaxis.grid(True, color=COLOR_GRID, linestyle="--", linewidth=0.6, zorder=0)
ax_b.set_axisbelow(True)

# Legend for Panel B using patch handles
cold_patch = mpatches.Patch(facecolor=COLOR_COLD, label="Cold cache")
warm_patch = mpatches.Patch(
    facecolor=COLOR_WARM, edgecolor=COLOR_COLD, linewidth=0.5, label="Warm cache"
)
ax_b.legend(handles=[cold_patch, warm_patch], loc="upper left",
            framealpha=0.85, edgecolor="#CCCCCC")

# Panel label (B)
ax_b.text(
    -0.14, 1.04, "(B)",
    transform=ax_b.transAxes,
    fontsize=PANEL_LABEL, fontweight="bold",
    va="top", ha="left",
)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------

base_path = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures/fig4_performance"

fig.savefig(base_path + ".svg", format="svg",
            bbox_inches="tight", dpi=300)
fig.savefig(base_path + ".pdf", format="pdf",
            bbox_inches="tight", dpi=300)

print(f"Saved: {base_path}.svg")
print(f"Saved: {base_path}.pdf")
