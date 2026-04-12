"""
Generate publication-quality system architecture diagram for ChemFuse.
Outputs: fig1_architecture.svg and fig1_architecture.pdf
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---------------------------------------------------------------------------
# Figure dimensions (180 mm x 140 mm at 300 dpi for journal)
# ---------------------------------------------------------------------------
MM_TO_INCH = 1 / 25.4
FIG_W = 180 * MM_TO_INCH   # ~7.09 in
FIG_H = 150 * MM_TO_INCH   # ~5.91 in (slightly taller for comfortable layout)

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.patch.set_facecolor("white")

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
COLOR_UI     = "#D6E9F8"   # light blue  – User Interface layer
COLOR_CORE   = "#D9F0D3"   # light green – Core Services layer
COLOR_COMP   = "#FCE8C8"   # light orange – Computation layer
COLOR_BORDER = "#4A4A4A"   # dark grey border
COLOR_ARROW  = "#2B6CB0"   # blue arrows (query downward)
COLOR_ARROW2 = "#276749"   # green arrows (results upward)
COLOR_SECTION_HEADER = "#2D3748"  # near-black for tier labels

# ---------------------------------------------------------------------------
# Font sizes
# ---------------------------------------------------------------------------
FS_TIER_LABEL  = 7.5   # tier header text
FS_BOX_MAIN    = 6.2   # primary label inside a box
FS_BOX_SUB     = 5.2   # secondary/subtitle text
FS_CAPTION     = 5.0   # footnote

# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------
TIER_X0 = 0.03          # left edge of all tier blocks
TIER_W  = 0.94          # width of tier blocks
TIER_RADIUS = 0.012     # rounded corner radius

# Tier Y positions and heights (bottom = 0, top = 1)
# Tiers from top to bottom: UI, Core, Computation
TIER_GAP  = 0.025
TIER_UI_Y0     = 0.74
TIER_UI_H      = 0.20
TIER_CORE_Y0   = 0.38
TIER_CORE_H    = 0.30
TIER_COMP_Y0   = 0.05
TIER_COMP_H    = 0.28

# ---------------------------------------------------------------------------
# Helper: draw a tier background rectangle with a left-side label bar
# ---------------------------------------------------------------------------
def draw_tier_bg(ax, x0, y0, w, h, color, label, text_color=COLOR_SECTION_HEADER):
    """Draw the tier background panel and rotated left-edge label."""
    # Main panel
    rect = FancyBboxPatch(
        (x0, y0), w, h,
        boxstyle=f"round,pad=0,rounding_size={TIER_RADIUS}",
        linewidth=0.8,
        edgecolor=COLOR_BORDER,
        facecolor=color,
        zorder=1,
    )
    ax.add_patch(rect)
    # Vertical label on left margin
    ax.text(
        x0 - 0.015, y0 + h / 2, label,
        ha="center", va="center",
        fontsize=FS_TIER_LABEL,
        fontweight="bold",
        color=text_color,
        rotation=90,
        zorder=5,
    )


# ---------------------------------------------------------------------------
# Helper: draw a box with optional two-line text
# ---------------------------------------------------------------------------
def draw_box(ax, x, y, w, h, label, sublabel=None,
             facecolor="white", edgecolor=COLOR_BORDER,
             lw=0.7, fontsize=FS_BOX_MAIN, sub_fontsize=FS_BOX_SUB,
             label_bold=False, zorder=3):
    """Draw a white rounded box with primary label and optional sublabel."""
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0,rounding_size=0.008",
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=zorder,
    )
    ax.add_patch(box)
    weight = "bold" if label_bold else "normal"
    if sublabel:
        ax.text(
            x + w / 2, y + h * 0.60, label,
            ha="center", va="center",
            fontsize=fontsize, fontweight=weight,
            color=COLOR_SECTION_HEADER, zorder=zorder + 1,
        )
        ax.text(
            x + w / 2, y + h * 0.26, sublabel,
            ha="center", va="center",
            fontsize=sub_fontsize,
            color="#555555", zorder=zorder + 1,
        )
    else:
        ax.text(
            x + w / 2, y + h / 2, label,
            ha="center", va="center",
            fontsize=fontsize, fontweight=weight,
            color=COLOR_SECTION_HEADER, zorder=zorder + 1,
        )


# ---------------------------------------------------------------------------
# Helper: draw a vertical arrow between tiers
# ---------------------------------------------------------------------------
def draw_arrow(ax, x, y_start, y_end, color, label="", lw=1.0):
    ax.annotate(
        "",
        xy=(x, y_end), xytext=(x, y_start),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            mutation_scale=7,
        ),
        zorder=6,
    )
    if label:
        mid_y = (y_start + y_end) / 2
        ax.text(
            x + 0.008, mid_y, label,
            ha="left", va="center",
            fontsize=4.5, color=color, zorder=7,
        )


# ===========================================================================
# TIER 1 – USER INTERFACE LAYER
# ===========================================================================
draw_tier_bg(ax, TIER_X0, TIER_UI_Y0, TIER_W, TIER_UI_H, COLOR_UI,
             "User Interface Layer")

# 5 boxes evenly across the tier
ui_items = [
    ("Python API",         "chemfuse.*"),
    ("R Package",          "(reticulate)"),
    ("CLI",                "(Click)"),
    ("Web UI",             "(Streamlit)"),
    ("Docker",             "container"),
]
N_UI = len(ui_items)
ui_box_w = 0.145
ui_box_h = 0.090
ui_y    = TIER_UI_Y0 + (TIER_UI_H - ui_box_h) / 2

# Compute x positions centred within the tier
ui_total_w = N_UI * ui_box_w + (N_UI - 1) * 0.015
ui_x_start = TIER_X0 + (TIER_W - ui_total_w) / 2
ui_box_xs = []

for i, (lbl, sub) in enumerate(ui_items):
    bx = ui_x_start + i * (ui_box_w + 0.015)
    ui_box_xs.append(bx)
    draw_box(ax, bx, ui_y, ui_box_w, ui_box_h,
             lbl, sub,
             facecolor="white",
             label_bold=True)


# ===========================================================================
# TIER 2 – CORE SERVICES
# ===========================================================================
draw_tier_bg(ax, TIER_X0, TIER_CORE_Y0, TIER_W, TIER_CORE_H, COLOR_CORE,
             "Core Services")

# ---- Row A: Database Adapters ----
db_label_y = TIER_CORE_Y0 + TIER_CORE_H - 0.040
ax.text(
    TIER_X0 + TIER_W / 2, db_label_y,
    "Database Adapters",
    ha="center", va="center",
    fontsize=FS_BOX_MAIN - 0.3,
    fontweight="bold",
    color=COLOR_SECTION_HEADER,
    zorder=5,
)

db_items = [
    "PubChem", "ChEMBL", "UniChem",
    "BindingDB", "Open Targets", "SureChEMBL*",
]
N_DB = len(db_items)
db_box_w = 0.128
db_box_h = 0.065
db_row_y  = TIER_CORE_Y0 + TIER_CORE_H - 0.118

db_total_w = N_DB * db_box_w + (N_DB - 1) * 0.010
db_x_start = TIER_X0 + (TIER_W - db_total_w) / 2
db_box_xs  = []

for i, name in enumerate(db_items):
    bx = db_x_start + i * (db_box_w + 0.010)
    db_box_xs.append(bx + db_box_w / 2)
    # SureChEMBL gets a slightly different edge to flag experimental
    is_exp = name.endswith("*")
    ec = "#CC7A00" if is_exp else COLOR_BORDER
    lw = 1.0 if is_exp else 0.7
    draw_box(ax, bx, db_row_y, db_box_w, db_box_h,
             name,
             facecolor="#FFFDF5" if is_exp else "white",
             edgecolor=ec, lw=lw,
             label_bold=False)

# ---- Row B: three service boxes ----
svc_items = [
    ("Async HTTP Client",    "httpx + asyncio"),
    ("SQLite Cache",         "TTL + LRU"),
    ("Export Engine",        "CSV · JSON · Excel · SDF"),
]
N_SVC = len(svc_items)
svc_box_w = 0.24
svc_box_h = 0.080
svc_row_y  = TIER_CORE_Y0 + 0.022

svc_total_w = N_SVC * svc_box_w + (N_SVC - 1) * 0.025
svc_x_start = TIER_X0 + (TIER_W - svc_total_w) / 2
svc_box_xs  = []

for i, (lbl, sub) in enumerate(svc_items):
    bx = svc_x_start + i * (svc_box_w + 0.025)
    svc_box_xs.append(bx + svc_box_w / 2)
    draw_box(ax, bx, svc_row_y, svc_box_w, svc_box_h,
             lbl, sub,
             facecolor="white", label_bold=True)


# ===========================================================================
# TIER 3 – COMPUTATION MODULE
# ===========================================================================
draw_tier_bg(ax, TIER_X0, TIER_COMP_Y0, TIER_W, TIER_COMP_H, COLOR_COMP,
             "Computation Module")

comp_items = [
    ("RDKit Descriptors &\nFingerprints",
     "200+ descriptors · 5 FP types"),
    ("Drug-likeness\nFilters",
     "Lipinski · Veber · Ghose\nEgan · Muegge · PAINS · QED"),
    ("ADMET Prediction",
     "admet-ai ML\n+ rule-based fallback"),
    ("Analysis",
     "Butina clustering\nUMAP / t-SNE / PCA · SAR"),
]
N_COMP = len(comp_items)
comp_box_w = 0.193
comp_box_h = 0.150
comp_y = TIER_COMP_Y0 + (TIER_COMP_H - comp_box_h) / 2 + 0.010

comp_total_w = N_COMP * comp_box_w + (N_COMP - 1) * 0.020
comp_x_start = TIER_X0 + (TIER_W - comp_total_w) / 2
comp_box_xs = []

for i, (lbl, sub) in enumerate(comp_items):
    bx = comp_x_start + i * (comp_box_w + 0.020)
    comp_box_xs.append(bx + comp_box_w / 2)
    draw_box(ax, bx, comp_y, comp_box_w, comp_box_h,
             lbl, sub,
             facecolor="white", label_bold=True,
             fontsize=FS_BOX_MAIN, sub_fontsize=FS_BOX_SUB)


# ===========================================================================
# DATA FLOW ARROWS  (between tiers)
# ===========================================================================

# --- UI -> Core (queries going down) ---
# Draw two arrows: one on left side, one on right side of the figure
# Left arrow: queries downward
arrow_lx = 0.155
arrow_rx = 0.845

# UI bottom -> Core top  (query direction)
draw_arrow(ax,
           arrow_lx,
           TIER_UI_Y0,                          # start (UI bottom)
           TIER_CORE_Y0 + TIER_CORE_H + 0.002,  # end   (Core top)
           color=COLOR_ARROW, label="queries", lw=1.1)

# Core top -> UI bottom  (results direction, offset slightly)
draw_arrow(ax,
           arrow_lx + 0.028,
           TIER_CORE_Y0 + TIER_CORE_H + 0.001,  # start
           TIER_UI_Y0,                            # end
           color=COLOR_ARROW2, label="results", lw=1.1)

# --- Core -> Computation (queries going down) ---
draw_arrow(ax,
           arrow_rx,
           TIER_CORE_Y0,                          # start (Core bottom)
           TIER_COMP_Y0 + TIER_COMP_H + 0.002,   # end   (Comp top)
           color=COLOR_ARROW, lw=1.1)

# Computation top -> Core bottom (results upward)
draw_arrow(ax,
           arrow_rx - 0.028,
           TIER_COMP_Y0 + TIER_COMP_H + 0.001,
           TIER_CORE_Y0,
           color=COLOR_ARROW2, lw=1.1)


# ===========================================================================
# LEGEND for arrows
# ===========================================================================
legend_x = 0.68
legend_y  = TIER_UI_Y0 + TIER_UI_H - 0.040

ax.annotate("", xy=(legend_x + 0.045, legend_y),
            xytext=(legend_x, legend_y),
            arrowprops=dict(arrowstyle="-|>", color=COLOR_ARROW,
                            lw=1.0, mutation_scale=6),
            zorder=8)
ax.text(legend_x + 0.050, legend_y, "query / request",
        va="center", ha="left", fontsize=4.8, color=COLOR_ARROW)

legend_y2 = legend_y - 0.025
ax.annotate("", xy=(legend_x + 0.045, legend_y2),
            xytext=(legend_x, legend_y2),
            arrowprops=dict(arrowstyle="-|>", color=COLOR_ARROW2,
                            lw=1.0, mutation_scale=6),
            zorder=8)
ax.text(legend_x + 0.050, legend_y2, "result / data",
        va="center", ha="left", fontsize=4.8, color=COLOR_ARROW2)


# ===========================================================================
# FIGURE TITLE
# ===========================================================================
ax.text(
    0.50, 0.975,
    "ChemFuse System Architecture",
    ha="center", va="center",
    fontsize=8.5, fontweight="bold",
    color=COLOR_SECTION_HEADER,
    transform=ax.transAxes,
    zorder=10,
)


# ===========================================================================
# FOOTNOTE
# ===========================================================================
ax.text(
    TIER_X0, TIER_COMP_Y0 - 0.022,
    "* SureChEMBL: experimental, limited functionality",
    ha="left", va="top",
    fontsize=FS_CAPTION,
    color="#666666",
    fontstyle="italic",
    zorder=10,
)


# ===========================================================================
# SAVE
# ===========================================================================
OUT_SVG = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures/fig1_architecture.svg"
OUT_PDF = "/home/juhur/PROJECTS/3INITIALSTAGE3/chemfuse/paper/figures/fig1_architecture.pdf"

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", dpi=300,
            facecolor="white", metadata={"Creator": "ChemFuse / matplotlib"})
fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight", dpi=300,
            facecolor="white")

print(f"Saved: {OUT_SVG}")
print(f"Saved: {OUT_PDF}")
plt.close(fig)
