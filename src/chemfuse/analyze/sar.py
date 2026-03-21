"""Structure-Activity Relationship (SAR) analysis.

Provides activity cliff detection and similarity network construction.
"""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.analyze.similarity import tanimoto_matrix
from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

nx = None
try:
    import networkx as nx  # type: ignore[import-not-found]

    _NX_AVAILABLE = True
except ImportError:
    _NX_AVAILABLE = False


def detect_activity_cliffs(
    compounds: list[dict[str, Any]],
    activity_col: str = "activity",
    sim_threshold: float = 0.8,
    fp_type: str = "morgan",
) -> list[dict[str, Any]]:
    """Detect activity cliffs in a set of compounds.

    An activity cliff is defined as a pair of compounds with high structural
    similarity but a large difference in activity.

    Args:
        compounds: List of dicts, each containing a 'smiles' key and the
            activity column specified by *activity_col*.
        activity_col: Name of the activity column in each compound dict.
        sim_threshold: Minimum Tanimoto similarity to consider a pair.
        fp_type: Fingerprint type for similarity computation.

    Returns:
        List of cliff dicts sorted by descending cliff_score. Each dict has:
        smiles_1, smiles_2, similarity, activity_1, activity_2,
        activity_diff, cliff_score.

    Raises:
        ValueError: If the activity column is missing from any compound.
    """
    if not compounds:
        return []

    # Validate activity column presence
    for i, c in enumerate(compounds):
        if activity_col not in c:
            raise ValueError(
                f"Compound at index {i} is missing activity column '{activity_col}'. "
                "Enrich compounds from ChEMBL or BindingDB first (SPEC-CF-003)."
            )

    smiles_list = [c["smiles"] for c in compounds]
    activities = [float(c[activity_col]) for c in compounds]

    sim_mat = tanimoto_matrix(smiles_list, fp_type=fp_type)

    cliffs: list[dict[str, Any]] = []
    n = len(compounds)
    for i in range(n):
        for j in range(i + 1, n):
            sim = float(sim_mat[i, j])
            if sim >= sim_threshold:
                act_diff = abs(activities[i] - activities[j])
                cliff_score = sim * act_diff
                cliffs.append({
                    "smiles_1": smiles_list[i],
                    "smiles_2": smiles_list[j],
                    "similarity": round(sim, 4),
                    "activity_1": activities[i],
                    "activity_2": activities[j],
                    "activity_diff": round(act_diff, 4),
                    "cliff_score": round(cliff_score, 4),
                })

    cliffs.sort(key=lambda c: c["cliff_score"], reverse=True)
    return cliffs


def build_similarity_network(
    compounds: list[dict[str, Any]],
    threshold: float = 0.5,
    fp_type: str = "morgan",
) -> object:
    """Build a networkx Graph of molecular similarity.

    Args:
        compounds: List of dicts, each containing a 'smiles' key.
        threshold: Minimum Tanimoto similarity for an edge.
        fp_type: Fingerprint type.

    Returns:
        networkx.Graph with nodes indexed 0..n-1 and edges for pairs
        exceeding *threshold*. Each node has 'smiles' and optional 'name'.

    Raises:
        OptionalDependencyError: If networkx is not installed.
    """
    if not _NX_AVAILABLE:
        raise OptionalDependencyError("networkx", "analyze")

    smiles_list = [c["smiles"] for c in compounds]
    sim_mat = tanimoto_matrix(smiles_list, fp_type=fp_type)

    G = nx.Graph()
    for i, c in enumerate(compounds):
        G.add_node(i, smiles=c["smiles"], name=c.get("name", c["smiles"]))

    n = len(smiles_list)
    for i in range(n):
        for j in range(i + 1, n):
            sim = float(sim_mat[i, j])
            if sim >= threshold:
                G.add_edge(i, j, weight=round(sim, 4))

    return G
