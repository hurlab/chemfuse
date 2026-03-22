"""
ChemFuse Fig. 2 - Multi-database Query Workflow Diagram
Publication-quality SVG/PDF for academic journal submission.
Horizontal (left-to-right) flowchart.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

# ---------------------------------------------------------------------------
# Canvas / figure setup
# ---------------------------------------------------------------------------
# 180 mm x 100 mm in inches (1 inch = 25.4 mm)
FIG_W = 180 / 25.4   # 7.087 inches
FIG_H = 100 / 25.4   # 3.937 inches

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 18)
ax.set_ylim(0, 10)
ax.axis("off")
fig.patch.set_facecolor("white")

# ---------------------------------------------------------------------------
# Color palette  (journal-friendly)
# ---------------------------------------------------------------------------
C_INPUT   = "#2C6FAC"   # deep blue  - input / output boxes
C_INPUT_L = "#D6E8F7"   # light blue fill
C_DB      = "#2A7D4F"   # green      - database query boxes
C_DB_L    = "#D4EDE1"   # light green fill
C_CACHE   = "#E8A020"   # amber      - cache diamond
C_CACHE_L = "#FEF3DC"   # light amber fill
C_COMP    = "#B84D00"   # burnt orange - computation
C_COMP_L  = "#FDEBD0"   # light orange fill
C_MERGE   = "#5C3F8F"   # purple     - merge
C_MERGE_L = "#EAE0F5"   # light purple fill
C_ARROW   = "#444444"
C_TEXT    = "#111111"

# ---------------------------------------------------------------------------
# Helper drawing functions
# ---------------------------------------------------------------------------

def rounded_box(ax, cx, cy, w, h, text, facecolor, edgecolor,
                fontsize=7.5, text_color=C_TEXT, bold=False,
                text_lines=None, zorder=3):
    """Draw a rounded rectangle centered at (cx, cy)."""
    x0, y0 = cx - w / 2, cy - h / 2
    box = FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle="round,pad=0.08",
        linewidth=0.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=zorder,
    )
    ax.add_patch(box)
    weight = "bold" if bold else "normal"
    if text_lines is None:
        ax.text(cx, cy, text, ha="center", va="center",
                fontsize=fontsize, color=text_color,
                fontweight=weight, zorder=zorder + 1,
                fontfamily="sans-serif")
    else:
        # Multi-line with explicit line spacing
        for i, line in enumerate(text_lines):
            offset = (len(text_lines) - 1) / 2 - i
            ax.text(cx, cy + offset * fontsize * 0.012,
                    line, ha="center", va="center",
                    fontsize=fontsize, color=text_color,
                    fontweight="bold" if (bold and i == 0) else weight,
                    zorder=zorder + 1, fontfamily="sans-serif")


def diamond(ax, cx, cy, w, h, text, facecolor, edgecolor,
            fontsize=6.5, zorder=3):
    """Draw a diamond shape centered at (cx, cy)."""
    pts = np.array([
        [cx,       cy + h / 2],
        [cx + w/2, cy        ],
        [cx,       cy - h / 2],
        [cx - w/2, cy        ],
    ])
    poly = plt.Polygon(pts, closed=True, linewidth=0.8,
                       edgecolor=edgecolor, facecolor=facecolor, zorder=zorder)
    ax.add_patch(poly)
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fontsize, color=C_TEXT, zorder=zorder + 1,
            fontfamily="sans-serif")


def arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=0.9, zorder=2,
          arrowstyle="-|>", mutation_scale=8, label=None, label_side="top"):
    """Draw a directional arrow."""
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=arrowstyle,
        color=color,
        linewidth=lw,
        mutation_scale=mutation_scale,
        zorder=zorder,
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arr)
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dy = 0.18 if label_side == "top" else -0.18
        ax.text(mx, my + dy, label, ha="center", va="center",
                fontsize=5.5, color="#555555", zorder=zorder + 1,
                fontfamily="sans-serif")


def bent_arrow(ax, points, color=C_ARROW, lw=0.9, zorder=2, arrowstyle="-|>"):
    """Draw a multi-segment arrow through a list of (x,y) waypoints."""
    from matplotlib.path import Path
    import matplotlib.patches as mpatches
    codes = [Path.MOVETO] + [Path.LINETO] * (len(points) - 1)
    path = Path(points, codes)
    patch = mpatches.FancyArrowPatch(
        path=path,
        arrowstyle=arrowstyle,
        color=color,
        linewidth=lw,
        mutation_scale=8,
        zorder=zorder,
    )
    ax.add_patch(patch)


def polyline(ax, points, color=C_ARROW, lw=0.9, zorder=2):
    """Draw a plain poly-line (no arrowhead) through waypoints."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.plot(xs, ys, color=color, lw=lw, zorder=zorder,
            solid_capstyle="round", solid_joinstyle="round")

# ---------------------------------------------------------------------------
# Layout constants  (x positions, y positions)
# ---------------------------------------------------------------------------
# X columns
X_IN      = 1.55   # Input query box center
X_GATHER  = 3.30   # asyncio.gather box center
X_CACHE   = 5.50   # Cache diamond center
X_DBBOX   = 7.50   # DB API box center
X_CACHEOK = 9.20   # "cache hit" small label / bypass join
X_MERGE   = 10.70  # Merge box center
X_ENR1    = 12.50  # Enrichment step 1
X_ENR2    = 13.90  # Enrichment step 2
X_ENR3    = 15.30  # Enrichment step 3
X_OUT     = 17.10  # Output box center

# Y rows for the three parallel database tracks
Y_TOP     = 7.40
Y_MID     = 5.00
Y_BOT     = 2.60
Y_GATHER  = (Y_TOP + Y_BOT) / 2   # 5.0 (same as mid for gather box)

# Enrichment Y  (placed at mid-level)
Y_ENR     = Y_MID

# Box dimensions
W_IN      = 2.70
H_IN      = 1.10
W_GATHER  = 1.30
H_GATHER  = 0.80
W_CACHE   = 1.20
H_CACHE   = 0.80
W_DB      = 1.70
H_DB      = 0.80
W_MERGE   = 1.60
H_MERGE   = 0.80
W_ENR     = 1.40
H_ENR     = 0.70
W_OUT     = 1.90
H_OUT     = 1.00

# Diamond half-widths/heights
DW = 0.85
DH = 0.60

# ---------------------------------------------------------------------------
# (1) Input box
# ---------------------------------------------------------------------------
IN_LINES = [
    "chemfuse.search(",
    "  'aspirin',",
    "  sources=['pubchem',",
    "   'chembl','bindingdb']",
    ")",
]
x0_in, y0_in = X_IN - W_IN / 2, Y_MID - 1.10
in_box = FancyBboxPatch(
    (x0_in, y0_in), W_IN, 2.20,
    boxstyle="round,pad=0.10",
    linewidth=1.0,
    edgecolor=C_INPUT,
    facecolor=C_INPUT_L,
    zorder=3,
)
ax.add_patch(in_box)
ax.text(X_IN, y0_in + 2.20 + 0.20, "User Query", ha="center", va="bottom",
        fontsize=7, fontweight="bold", color=C_INPUT,
        fontfamily="sans-serif", zorder=4)
for i, line in enumerate(IN_LINES):
    ax.text(X_IN, y0_in + 1.80 - i * 0.38, line,
            ha="center", va="center",
            fontsize=6.2, color="#1A1A1A", family="monospace", zorder=4)

# ---------------------------------------------------------------------------
# (2) asyncio.gather() node
# ---------------------------------------------------------------------------
rounded_box(ax, X_GATHER, Y_MID, W_GATHER, H_GATHER,
            "asyncio\n.gather()",
            C_DB_L, C_DB, fontsize=6.5, bold=True)

# Arrow: input → gather
arrow(ax, X_IN + W_IN / 2, Y_MID,
          X_GATHER - W_GATHER / 2, Y_MID)

# Fan-out arrows from gather to 3 tracks
gather_right = X_GATHER + W_GATHER / 2
FAN_X = 3.95   # intermediate x for fan-out elbow

# top track
polyline(ax, [(gather_right, Y_MID),
              (FAN_X, Y_MID),
              (FAN_X, Y_TOP)], color=C_ARROW)
arrow(ax, FAN_X, Y_TOP, X_CACHE - DW, Y_TOP)

# mid track
arrow(ax, gather_right, Y_MID, X_CACHE - DW, Y_MID)

# bottom track
polyline(ax, [(gather_right, Y_MID),
              (FAN_X, Y_MID),
              (FAN_X, Y_BOT)], color=C_ARROW)
arrow(ax, FAN_X, Y_BOT, X_CACHE - DW, Y_BOT)

# ---------------------------------------------------------------------------
# (3) Cache diamonds + DB boxes for each track
# ---------------------------------------------------------------------------
tracks = [
    ("PubChem\nPUG-REST", "compounds", Y_TOP),
    ("ChEMBL\nREST API",  "compounds", Y_MID),
    ("BindingDB\nSOAP",   "binding data", Y_BOT),
]

CACHE_JOIN_X = 9.50   # x where "yes" and "no" paths rejoin

for db_label, return_label, yk in tracks:
    # Cache diamond
    diamond(ax, X_CACHE, yk, DW * 2, DH * 2, "Cached?",
            C_CACHE_L, C_CACHE, fontsize=6.5)

    # DB box (No path → right of diamond)
    rounded_box(ax, X_DBBOX, yk, W_DB, H_DB,
                db_label, C_DB_L, C_DB, fontsize=6.8, bold=True)

    # No arrow: diamond right → DB box
    arrow(ax, X_CACHE + DW, yk, X_DBBOX - W_DB / 2, yk,
          label="No", label_side="top")

    # DB box → merge junction
    arrow(ax, X_DBBOX + W_DB / 2, yk, CACHE_JOIN_X, yk,
          label=f"→ {return_label}", label_side="top")

    # Yes path: diamond top → arc over → join x
    # route: up from diamond top, bend right to join
    YES_BYPASS_Y = yk + DH + 0.55
    polyline(ax, [(X_CACHE, yk + DH),
                  (X_CACHE, YES_BYPASS_Y),
                  (CACHE_JOIN_X, YES_BYPASS_Y)], color=C_CACHE, lw=0.85)
    arrow(ax, CACHE_JOIN_X, YES_BYPASS_Y, CACHE_JOIN_X, yk + 0.12,
          color=C_CACHE, label="Yes", label_side="top")

    ax.text(X_CACHE - DW - 0.08, yk + 0.20, "cache\nhit",
            ha="right", va="center",
            fontsize=5.5, color=C_CACHE, fontfamily="sans-serif")

# ---------------------------------------------------------------------------
# (4) Merge box
# ---------------------------------------------------------------------------
rounded_box(ax, X_MERGE, Y_MID, W_MERGE, H_MERGE,
            "InChIKey-based\nMerge",
            C_MERGE_L, C_MERGE, fontsize=6.8, bold=True)

# Fan-in arrows from three CACHE_JOIN_X points to merge box
MERGE_ELBOW_X = 10.05
for yk in [Y_TOP, Y_MID, Y_BOT]:
    if yk == Y_MID:
        arrow(ax, CACHE_JOIN_X, yk, X_MERGE - W_MERGE / 2, yk)
    else:
        polyline(ax, [(CACHE_JOIN_X, yk),
                      (MERGE_ELBOW_X, yk),
                      (MERGE_ELBOW_X, Y_MID)], color=C_ARROW)
arrow(ax, MERGE_ELBOW_X, Y_MID, X_MERGE - W_MERGE / 2, Y_MID)

# ---------------------------------------------------------------------------
# (5) Enrichment pipeline (optional, shown as sequential boxes)
# ---------------------------------------------------------------------------
enrichment = [
    (X_ENR1, "Descriptors\n200+ RDKit"),
    (X_ENR2, "Drug-likeness\n7 filters"),
    (X_ENR3, "ADMET\nprediction"),
]

# Arrow from merge to first enrichment
arrow(ax, X_MERGE + W_MERGE / 2, Y_ENR,
          X_ENR1 - W_ENR / 2, Y_ENR,
          label="enrich", label_side="top")

for i, (xc, label) in enumerate(enrichment):
    rounded_box(ax, xc, Y_ENR, W_ENR, H_ENR,
                label, C_COMP_L, C_COMP, fontsize=6.5, bold=True)
    if i < len(enrichment) - 1:
        arrow(ax, xc + W_ENR / 2, Y_ENR,
                  enrichment[i + 1][0] - W_ENR / 2, Y_ENR)

# Optional bracket
BRK_Y0 = Y_ENR - H_ENR / 2 - 0.30
BRK_Y1 = Y_ENR + H_ENR / 2 + 0.30
BRK_X0 = X_ENR1 - W_ENR / 2 - 0.10
BRK_X1 = X_ENR3 + W_ENR / 2 + 0.10
for bx in [BRK_X0, BRK_X1]:
    polyline(ax, [(bx, BRK_Y0), (bx, BRK_Y1)],
             color="#888888", lw=0.6)
polyline(ax, [(BRK_X0, BRK_Y0), (BRK_X1, BRK_Y0)],
         color="#888888", lw=0.6)
polyline(ax, [(BRK_X0, BRK_Y1), (BRK_X1, BRK_Y1)],
         color="#888888", lw=0.6)
ax.text((BRK_X0 + BRK_X1) / 2, BRK_Y0 - 0.22,
        "optional enrichment pipeline",
        ha="center", va="top", fontsize=5.5,
        color="#888888", fontfamily="sans-serif", style="italic")

# ---------------------------------------------------------------------------
# (6) Output box (CompoundCollection)
# ---------------------------------------------------------------------------
OUT_LINES = [
    "CompoundCollection",
    "━━━━━━━━━━━━━━━━",
    "unified data",
]
out_x0 = X_OUT - W_OUT / 2
out_y0 = Y_MID - H_OUT / 2 - 0.10
out_box = FancyBboxPatch(
    (out_x0, out_y0), W_OUT, H_OUT + 0.20,
    boxstyle="round,pad=0.10",
    linewidth=1.0,
    edgecolor=C_INPUT,
    facecolor=C_INPUT_L,
    zorder=3,
)
ax.add_patch(out_box)
ax.text(X_OUT, out_y0 + H_OUT + 0.20 + 0.20,
        "Output", ha="center", va="bottom",
        fontsize=7, fontweight="bold", color=C_INPUT,
        fontfamily="sans-serif", zorder=4)
for i, line in enumerate(OUT_LINES):
    fw = "bold" if i == 0 else "normal"
    ax.text(X_OUT, out_y0 + H_OUT + 0.20 - 0.22 - i * 0.32,
            line, ha="center", va="center",
            fontsize=6.5 if i == 0 else 6.0,
            fontweight=fw, color="#1A1A1A",
            fontfamily="sans-serif" if i != 1 else "monospace",
            zorder=4)

# Arrow from last enrichment to output
arrow(ax, X_ENR3 + W_ENR / 2, Y_ENR,
          out_x0, Y_MID)

# ---------------------------------------------------------------------------
# (7) Section labels (column headers)
# ---------------------------------------------------------------------------
headers = [
    (X_IN,     "Input"),
    (X_GATHER, "Dispatch"),
    ((X_CACHE + X_DBBOX) / 2, "Cache / Fetch"),
    (X_MERGE,  "Merge"),
    ((X_ENR1 + X_ENR3) / 2, "Enrichment"),
    (X_OUT,    "Output"),
]
HEADER_Y = 9.55
for hx, ht in headers:
    ax.text(hx, HEADER_Y, ht,
            ha="center", va="center",
            fontsize=6.5, color="#555555",
            fontfamily="sans-serif", fontstyle="italic",
            zorder=4)
# Thin separator line
ax.axhline(HEADER_Y - 0.30, color="#CCCCCC", lw=0.5, zorder=1)

# ---------------------------------------------------------------------------
# (8) Legend
# ---------------------------------------------------------------------------
legend_items = [
    (C_INPUT_L, C_INPUT, "Input / Output"),
    (C_DB_L,    C_DB,    "Database Query"),
    (C_CACHE_L, C_CACHE, "Cache Check"),
    (C_COMP_L,  C_COMP,  "Computation"),
    (C_MERGE_L, C_MERGE, "Data Merge"),
]
LEG_X = 0.15
LEG_Y = 1.60
ax.text(LEG_X, LEG_Y + 0.55, "Legend",
        fontsize=6.5, fontweight="bold",
        color="#333333", fontfamily="sans-serif", zorder=4)
for i, (fc, ec, label) in enumerate(legend_items):
    bx = LEG_X
    by = LEG_Y - i * 0.40
    r = FancyBboxPatch((bx, by - 0.12), 0.35, 0.24,
                       boxstyle="round,pad=0.04",
                       linewidth=0.7, edgecolor=ec, facecolor=fc, zorder=4)
    ax.add_patch(r)
    ax.text(bx + 0.45, by, label,
            fontsize=6.0, va="center",
            color="#333333", fontfamily="sans-serif", zorder=4)

# ---------------------------------------------------------------------------
# (9) Figure title / caption line
# ---------------------------------------------------------------------------
ax.text(9.0, 0.28,
        "Fig. 2  |  ChemFuse multi-database query workflow",
        ha="center", va="center",
        fontsize=7.5, color="#222222",
        fontfamily="sans-serif", fontweight="bold", zorder=4)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
OUT_SVG = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures/fig2_workflow.svg"
OUT_PDF = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures/fig2_workflow.pdf"

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight",
            dpi=300, facecolor="white")
fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight",
            dpi=300, facecolor="white")

print(f"Saved SVG: {OUT_SVG}")
print(f"Saved PDF: {OUT_PDF}")
plt.close(fig)
