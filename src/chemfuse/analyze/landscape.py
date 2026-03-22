"""SAR landscape visualization: activity-colored chemical space with cliff edges."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

# Lazy imports — kept at module scope so tests can patch them easily.
# Both are guarded by OptionalDependencyError inside their own modules.
from chemfuse.analyze.chemspace import reduce_dimensions
from chemfuse.analyze.similarity import tanimoto_matrix

logger = logging.getLogger(__name__)


def sar_landscape(
    smiles_list: list[str],
    activities: list[float],
    method: str = "umap",
    fp_type: str = "morgan",
    cliff_threshold: float = 0.8,
    activity_diff_threshold: float = 1.0,
    colormap: str = "RdYlGn_r",
) -> dict[str, Any]:
    """Generate SAR landscape data with chemical space coordinates and activity cliff edges.

    Args:
        smiles_list: SMILES strings for compounds.
        activities: Activity values (e.g., pIC50) for each compound.
        method: Dimensionality reduction method ("umap", "tsne", "pca").
        fp_type: Fingerprint type for similarity.
        cliff_threshold: Minimum Tanimoto similarity for cliff detection.
        activity_diff_threshold: Minimum activity difference for cliff detection.
        colormap: Matplotlib/Plotly colormap name.

    Returns:
        dict with keys:
            "coords": np.ndarray shape (n, 2) -- 2D coordinates
            "activities": list[float] -- activity values (same order)
            "smiles": list[str] -- SMILES (same order)
            "cliff_pairs": list[dict] -- each dict has "i", "j", "similarity", "activity_diff"
            "figure": plotly Figure object (if plotly available, else None)

    Raises:
        ValueError: If smiles_list and activities have different lengths.
        ValueError: If smiles_list is empty.
    """
    if len(smiles_list) != len(activities):
        raise ValueError(
            f"smiles_list (len={len(smiles_list)}) and activities "
            f"(len={len(activities)}) must have the same length."
        )

    if not smiles_list:
        return {
            "coords": np.empty((0, 2), dtype=float),
            "activities": [],
            "smiles": [],
            "cliff_pairs": [],
            "figure": None,
        }

    # 1. Reduce dimensions to 2D coordinates
    coords = reduce_dimensions(smiles_list, method=method, fp_type=fp_type, n_components=2)

    # 2. Detect activity cliffs
    cliff_pairs = _detect_cliff_pairs(
        smiles_list,
        activities,
        fp_type=fp_type,
        cliff_threshold=cliff_threshold,
        activity_diff_threshold=activity_diff_threshold,
    )

    # 3. Build optional Plotly figure
    figure = _build_plotly_figure(
        coords=coords,
        activities=activities,
        smiles_list=smiles_list,
        cliff_pairs=cliff_pairs,
        colormap=colormap,
        method=method,
    )

    return {
        "coords": coords,
        "activities": list(activities),
        "smiles": list(smiles_list),
        "cliff_pairs": cliff_pairs,
        "figure": figure,
    }


def _detect_cliff_pairs(
    smiles_list: list[str],
    activities: list[float],
    fp_type: str,
    cliff_threshold: float,
    activity_diff_threshold: float,
) -> list[dict[str, Any]]:
    """Detect activity cliff pairs using Tanimoto similarity and activity difference.

    Args:
        smiles_list: SMILES strings.
        activities: Activity values aligned with smiles_list.
        fp_type: Fingerprint type for similarity computation.
        cliff_threshold: Minimum Tanimoto similarity to qualify as a cliff pair.
        activity_diff_threshold: Minimum activity difference to qualify as a cliff.

    Returns:
        List of cliff pair dicts with keys "i", "j", "similarity", "activity_diff".
    """
    n = len(smiles_list)
    if n < 2:
        return []

    sim_mat = tanimoto_matrix(smiles_list, fp_type=fp_type)
    cliff_pairs: list[dict[str, Any]] = []

    for i in range(n):
        for j in range(i + 1, n):
            sim = float(sim_mat[i, j])
            act_diff = abs(float(activities[i]) - float(activities[j]))
            if sim >= cliff_threshold and act_diff >= activity_diff_threshold:
                cliff_pairs.append({
                    "i": i,
                    "j": j,
                    "similarity": round(sim, 4),
                    "activity_diff": round(act_diff, 4),
                })

    return cliff_pairs


def _build_plotly_figure(
    coords: np.ndarray,
    activities: list[float],
    smiles_list: list[str],
    cliff_pairs: list[dict[str, Any]],
    colormap: str,
    method: str,
) -> Any:
    """Build a Plotly scatter plot with activity coloring and cliff edge lines.

    Args:
        coords: 2D coordinates array of shape (n, 2).
        activities: Activity values for coloring.
        smiles_list: SMILES strings for hover text.
        cliff_pairs: List of cliff pair dicts with "i" and "j" indices.
        colormap: Plotly colorscale name.
        method: Dimensionality reduction method name (used for axis labels).

    Returns:
        plotly.graph_objects.Figure, or None if plotly is not installed.
    """
    try:
        import plotly.graph_objects as go  # type: ignore[import-not-found]
    except ImportError:
        logger.debug("plotly not installed; figure will be None.")
        return None

    fig = go.Figure()

    # Scatter points colored by activity
    fig.add_trace(go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode="markers",
        marker=dict(
            color=activities,
            colorscale=colormap,
            showscale=True,
            colorbar=dict(title="Activity"),
            size=8,
        ),
        text=[
            f"SMILES: {s}<br>Activity: {a:.2f}"
            for s, a in zip(smiles_list, activities, strict=True)
        ],
        hoverinfo="text",
        name="Compounds",
    ))

    # Cliff edges as dotted red lines
    for pair in cliff_pairs:
        i, j = pair["i"], pair["j"]
        fig.add_trace(go.Scatter(
            x=[coords[i, 0], coords[j, 0]],
            y=[coords[i, 1], coords[j, 1]],
            mode="lines",
            line=dict(color="red", width=1, dash="dot"),
            showlegend=False,
            hoverinfo="skip",
        ))

    method_label = method.upper()
    fig.update_layout(
        title="SAR Landscape",
        xaxis_title=f"{method_label} 1",
        yaxis_title=f"{method_label} 2",
        template="plotly_white",
    )
    return fig
