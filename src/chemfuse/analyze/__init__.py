"""ChemFuse analysis tools: clustering, chemical space, SAR, similarity, scaffolds, and diversity."""

from __future__ import annotations

from chemfuse.analyze.chemspace import plot_chemical_space, reduce_dimensions
from chemfuse.analyze.clustering import (
    butina_clustering,
    cluster_compounds,
    compute_silhouette,
    kmeans_clustering,
)
from chemfuse.analyze.diversity import diversity_score, maxmin_pick
from chemfuse.analyze.sar import build_similarity_network, detect_activity_cliffs
from chemfuse.analyze.scaffolds import (
    generic_scaffold,
    group_by_scaffold,
    murcko_scaffold,
    scaffold_frequency,
)
from chemfuse.analyze.similarity import (
    bulk_tanimoto,
    substructure_match_atoms,
    substructure_search,
    tanimoto_matrix,
    tanimoto_similarity,
)

__all__ = [
    "butina_clustering",
    "cluster_compounds",
    "compute_silhouette",
    "kmeans_clustering",
    "reduce_dimensions",
    "plot_chemical_space",
    "detect_activity_cliffs",
    "build_similarity_network",
    "tanimoto_similarity",
    "tanimoto_matrix",
    "bulk_tanimoto",
    "murcko_scaffold",
    "generic_scaffold",
    "scaffold_frequency",
    "group_by_scaffold",
    "substructure_search",
    "substructure_match_atoms",
    "maxmin_pick",
    "diversity_score",
]
