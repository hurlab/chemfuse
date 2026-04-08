"""Molecular clustering: Butina (Taylor-Butina) and KMeans algorithms.

Requires rdkit (for Butina) and scikit-learn (for KMeans).
"""

from __future__ import annotations

import logging

import numpy as np

from chemfuse.compute.fingerprints import fp_matrix as _shared_fp_matrix
from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

Chem = None
DataStructs = None
MACCSkeys = None
Butina = None
try:
    from rdkit import Chem, DataStructs  # type: ignore[import-not-found]
    from rdkit.Chem import MACCSkeys, rdFingerprintGenerator  # type: ignore[import-not-found]
    from rdkit.ML.Cluster import Butina  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
    _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
except ImportError:
    _RDKIT_AVAILABLE = False

KMeans = None
silhouette_score = None
try:
    from sklearn.cluster import KMeans  # type: ignore[import-not-found]
    from sklearn.metrics import silhouette_score  # type: ignore[import-not-found]

    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False


def _require_rdkit() -> None:
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "compute")


def _require_sklearn() -> None:
    if not _SKLEARN_AVAILABLE:
        raise OptionalDependencyError("scikit-learn", "analyze")


def _get_fp(smiles: str, fp_type: str = "morgan") -> object:
    """Return an RDKit fingerprint object for the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    if fp_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    if fp_type == "rdkit":
        return Chem.RDKFingerprint(mol)
    return _MORGAN_GEN.GetFingerprint(mol)


def _fp_matrix(smiles_list: list[str], fp_type: str = "morgan", n_bits: int = 2048) -> np.ndarray:
    """Build a numpy bit-vector matrix from SMILES.

    Delegates to the shared fp_matrix() in chemfuse.compute.fingerprints
    which correctly handles MACCS (167-bit) and other fingerprint dimensions.
    """
    return _shared_fp_matrix(smiles_list, fp_type=fp_type, n_bits=n_bits)


def butina_clustering(
    smiles_list: list[str],
    cutoff: float = 0.4,
    fp_type: str = "morgan",
) -> list[int]:
    """Cluster molecules using the Taylor-Butina algorithm.

    Args:
        smiles_list: List of SMILES strings.
        cutoff: Distance cutoff (1 - Tanimoto). Smaller = tighter clusters.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

    Returns:
        List of integer cluster labels, one per input SMILES.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles_list:
        return []

    fps = []
    valid_indices: list[int] = []
    for i, smi in enumerate(smiles_list):
        try:
            fps.append(_get_fp(smi, fp_type))
            valid_indices.append(i)
        except (ValueError, Exception) as exc:
            logger.warning("butina_clustering: skipping %r: %s", smi, exc)

    n = len(fps)
    if n < 2:
        return list(range(len(smiles_list)))

    # Condensed distance matrix for Butina
    dists: list[float] = []
    for i in range(1, n):
        for j in range(i):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            dists.append(1.0 - sim)

    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)

    labels = [-1] * len(smiles_list)
    for cluster_id, cluster in enumerate(clusters):
        for member_idx in cluster:
            labels[valid_indices[member_idx]] = cluster_id

    # Assign unprocessed (invalid SMILES) to their own clusters
    next_label = max(labels) + 1 if any(x >= 0 for x in labels) else 0
    for i in range(len(labels)):
        if labels[i] == -1:
            labels[i] = next_label
            next_label += 1

    return labels


def kmeans_clustering(
    smiles_list: list[str],
    n_clusters: int = 5,
    fp_type: str = "morgan",
) -> list[int]:
    """Cluster molecules using KMeans on fingerprint bit-vectors.

    Args:
        smiles_list: List of SMILES strings.
        n_clusters: Number of clusters.
        fp_type: Fingerprint type.

    Returns:
        List of integer cluster labels, one per input SMILES.

    Raises:
        OptionalDependencyError: If scikit-learn is not installed.
    """
    _require_rdkit()
    _require_sklearn()

    if not smiles_list:
        return []

    fp_mat = _fp_matrix(smiles_list, fp_type=fp_type)
    actual_k = min(n_clusters, len(smiles_list))
    if actual_k < 2:
        return [0] * len(smiles_list)

    km = KMeans(n_clusters=actual_k, random_state=42, n_init=10)
    return km.fit_predict(fp_mat).tolist()


def compute_silhouette(
    smiles_list: list[str],
    labels: list[int],
    fp_type: str = "morgan",
) -> float | None:
    """Compute the silhouette score for a clustering result.

    Args:
        smiles_list: List of SMILES strings.
        labels: Cluster labels (one per SMILES).
        fp_type: Fingerprint type.

    Returns:
        Silhouette score in [-1, 1], or None if scikit-learn is unavailable
        or only one cluster exists.
    """
    if not _SKLEARN_AVAILABLE:
        return None
    if not _RDKIT_AVAILABLE:
        return None
    n_unique = len(set(labels))
    if n_unique < 2 or n_unique >= len(smiles_list):
        return None
    fp_mat = _fp_matrix(smiles_list, fp_type=fp_type)
    try:
        return float(silhouette_score(fp_mat, labels))
    except Exception as exc:
        logger.warning("silhouette_score failed: %s", exc)
        return None


def cluster_compounds(
    smiles_list: list[str],
    method: str = "butina",
    n_clusters: int = 5,
    threshold: float = 0.65,
    fp_type: str = "morgan",
    n_bits: int = 2048,
) -> list[int]:
    """Cluster compounds by structural similarity. Convenience wrapper.

    Args:
        smiles_list: List of SMILES strings.
        method: Clustering method ('butina' or 'kmeans').
        n_clusters: Number of clusters for KMeans.
        threshold: Distance cutoff for Butina (1 - Tanimoto similarity).
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').
        n_bits: Fingerprint bit length.

    Returns:
        List of integer cluster labels, one per input SMILES.
    """
    if method == "kmeans":
        return kmeans_clustering(smiles_list, n_clusters=n_clusters, fp_type=fp_type)
    return butina_clustering(smiles_list, cutoff=1.0 - threshold, fp_type=fp_type)


__all__ = [
    "butina_clustering",
    "cluster_compounds",
    "compute_silhouette",
    "kmeans_clustering",
]
