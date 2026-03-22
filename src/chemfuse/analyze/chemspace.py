"""Chemical space dimensionality reduction and visualization.

Supports UMAP (preferred), t-SNE, and PCA for 2D projection of fingerprint
feature vectors. Requires scikit-learn; umap-learn is optional (falls back
to PCA when absent).
"""

from __future__ import annotations

import logging
import warnings
from typing import Any

import numpy as np

from chemfuse.compute.fingerprints import RDKIT_AVAILABLE as _RDKIT_AVAILABLE
from chemfuse.compute.fingerprints import fp_matrix as _shared_fp_matrix
from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

PCA = None
TSNE = None
try:
    from sklearn.decomposition import PCA  # type: ignore[import-not-found]
    from sklearn.manifold import TSNE  # type: ignore[import-not-found]

    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False

UMAP = None
try:
    from umap import UMAP  # type: ignore[import-not-found]

    _UMAP_AVAILABLE = True
except ImportError:
    _UMAP_AVAILABLE = False
    logger.debug("umap-learn not installed; UMAP will fall back to PCA.")


def _require_rdkit() -> None:
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "compute")


def _require_sklearn() -> None:
    if not _SKLEARN_AVAILABLE:
        raise OptionalDependencyError("scikit-learn", "analyze")


def _fp_matrix(smiles_list: list[str], fp_type: str = "morgan", n_bits: int = 2048) -> np.ndarray:
    """Build a numpy bit-vector matrix from SMILES.

    Delegates to the shared fp_matrix() in chemfuse.compute.fingerprints
    which correctly handles MACCS (167-bit) and other fingerprint dimensions.
    """
    return _shared_fp_matrix(smiles_list, fp_type=fp_type, n_bits=n_bits)


def reduce_dimensions(
    smiles_list: list[str],
    method: str = "umap",
    fp_type: str = "morgan",
    n_bits: int = 2048,
    n_components: int = 2,
    random_state: int = 42,
    **kwargs: Any,
) -> np.ndarray:
    """Reduce molecular fingerprints to 2D (or nD) coordinates.

    Args:
        smiles_list: List of SMILES strings.
        method: Reduction method: 'umap', 'tsne', or 'pca'.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').
        n_bits: Fingerprint bit length (for morgan/rdkit).
        n_components: Number of output dimensions (default 2).
        random_state: Random seed for reproducibility.
        **kwargs: Additional arguments passed to the underlying reducer.

    Returns:
        numpy array of shape (n_compounds, n_components).

    Raises:
        OptionalDependencyError: If required packages are missing.
    """
    _require_rdkit()
    _require_sklearn()

    fp_mat = _fp_matrix(smiles_list, fp_type=fp_type, n_bits=n_bits)
    n_samples = fp_mat.shape[0]

    if n_samples < 2:
        return np.zeros((n_samples, n_components))

    method = method.lower()

    if method == "umap":
        if not _UMAP_AVAILABLE:
            warnings.warn(
                "umap-learn is not installed; falling back to PCA. "
                "Install with: pip install umap-learn",
                UserWarning,
                stacklevel=2,
            )
            method = "pca"
        else:
            n_neighbors = min(kwargs.get("n_neighbors", 15), n_samples - 1)
            reducer = UMAP(
                n_components=n_components,
                n_neighbors=max(2, n_neighbors),
                min_dist=kwargs.get("min_dist", 0.1),
                random_state=random_state,
                metric="jaccard",
            )
            return reducer.fit_transform(fp_mat)

    if method == "tsne":
        perplexity = min(kwargs.get("perplexity", 30), max(1, (n_samples - 1) // 3))
        reducer = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            random_state=random_state,
            max_iter=1000,
        )
        return reducer.fit_transform(fp_mat)

    # PCA (default / fallback)
    eff_components = min(n_components, n_samples, fp_mat.shape[1])
    reducer_pca = PCA(n_components=eff_components, random_state=random_state)
    coords = reducer_pca.fit_transform(fp_mat)
    if coords.shape[1] < n_components:
        padding = np.zeros((coords.shape[0], n_components - coords.shape[1]))
        coords = np.column_stack([coords, padding])
    return coords


def plot_chemical_space(
    coordinates: np.ndarray,
    labels: list[str] | None = None,
    colors: list[Any] | None = None,
    title: str = "Chemical Space",
) -> object:
    """Create an interactive plotly Figure for 2D chemical space.

    Args:
        coordinates: Array of shape (n, 2) with 2D coordinates.
        labels: Optional hover labels.
        colors: Optional values or strings for point coloring.
        title: Plot title.

    Returns:
        plotly.graph_objects.Figure.

    Raises:
        OptionalDependencyError: If plotly is not installed.
    """
    try:
        import plotly.graph_objects as go  # type: ignore[import-not-found]
    except ImportError as exc:
        raise OptionalDependencyError("plotly", "viz") from exc

    n = len(coordinates)
    hover_text = labels if labels else [f"Compound {i}" for i in range(n)]

    marker: dict[str, Any] = {"size": 8}
    if colors is not None:
        if all(isinstance(c, str) for c in colors):
            marker["color"] = colors
        else:
            marker["color"] = colors
            marker["colorscale"] = "Viridis"
            marker["showscale"] = True

    fig = go.Figure(
        data=go.Scatter(
            x=coordinates[:, 0],
            y=coordinates[:, 1],
            mode="markers",
            marker=marker,
            text=hover_text,
            hoverinfo="text",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Dim 1",
        yaxis_title="Dim 2",
        height=500,
    )
    return fig
