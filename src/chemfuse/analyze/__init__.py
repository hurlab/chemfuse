"""ChemFuse analysis tools: clustering, chemical space, SAR, and similarity."""

from __future__ import annotations

from chemfuse.analyze.chemspace import plot_chemical_space, reduce_dimensions
from chemfuse.analyze.clustering import butina_clustering, kmeans_clustering
from chemfuse.analyze.sar import build_similarity_network, detect_activity_cliffs
from chemfuse.analyze.similarity import bulk_tanimoto, tanimoto_matrix, tanimoto_similarity

__all__ = [
    "butina_clustering",
    "kmeans_clustering",
    "reduce_dimensions",
    "plot_chemical_space",
    "detect_activity_cliffs",
    "build_similarity_network",
    "tanimoto_similarity",
    "tanimoto_matrix",
    "bulk_tanimoto",
]
