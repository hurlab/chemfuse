"""Local Tanimoto similarity computation using RDKit fingerprints.

No external API calls - all computation is local.
"""

from __future__ import annotations

import logging

import numpy as np

from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

Chem = None
DataStructs = None
AllChem = None
MACCSkeys = None
try:
    from rdkit import Chem, DataStructs  # type: ignore[import-not-found]
    from rdkit.Chem import AllChem, MACCSkeys  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


def _require_rdkit() -> None:
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "compute")


def _get_fingerprint(smiles: str, fp_type: str = "morgan") -> object:
    """Compute an RDKit fingerprint for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    if fp_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    if fp_type == "rdkit":
        return Chem.RDKFingerprint(mol)
    # Default: Morgan (radius=2, 2048 bits)
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def tanimoto_similarity(
    smiles_a: str,
    smiles_b: str,
    fp_type: str = "morgan",
) -> float:
    """Compute the Tanimoto similarity between two SMILES strings.

    Args:
        smiles_a: First SMILES string.
        smiles_b: Second SMILES string.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

    Returns:
        Tanimoto similarity in [0, 1].

    Raises:
        OptionalDependencyError: If RDKit is not installed.
        ValueError: If either SMILES is invalid.
    """
    _require_rdkit()
    fp_a = _get_fingerprint(smiles_a, fp_type)
    fp_b = _get_fingerprint(smiles_b, fp_type)
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def bulk_tanimoto(
    query: str,
    targets: list[str],
    fp_type: str = "morgan",
) -> list[tuple[str, float]]:
    """Compute Tanimoto similarity between a query and multiple target SMILES.

    Args:
        query: Query SMILES string.
        targets: List of target SMILES strings.
        fp_type: Fingerprint type.

    Returns:
        List of (smiles, similarity) tuples sorted by descending similarity.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()
    qfp = _get_fingerprint(query, fp_type)
    results: list[tuple[str, float]] = []
    for smi in targets:
        try:
            tfp = _get_fingerprint(smi, fp_type)
            sim = float(DataStructs.TanimotoSimilarity(qfp, tfp))
        except (ValueError, Exception) as exc:
            logger.warning("bulk_tanimoto: skipping %r: %s", smi, exc)
            sim = float("nan")
        results.append((smi, sim))
    results.sort(key=lambda x: x[1], reverse=True)
    return results


def tanimoto_matrix(
    smiles_list: list[str],
    fp_type: str = "morgan",
) -> np.ndarray:
    """Compute the full pairwise Tanimoto similarity matrix.

    Args:
        smiles_list: List of SMILES strings.
        fp_type: Fingerprint type.

    Returns:
        2D numpy array of shape (n, n) with pairwise similarities.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()
    n = len(smiles_list)
    fps = []
    for smi in smiles_list:
        try:
            fps.append(_get_fingerprint(smi, fp_type))
        except (ValueError, Exception) as exc:
            logger.warning("tanimoto_matrix: could not compute FP for %r: %s", smi, exc)
            fps.append(None)

    matrix = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if fps[i] is None or fps[j] is None:
                sim = float("nan")
            else:
                sim = float(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            matrix[i, j] = sim
            matrix[j, i] = sim
    return matrix


def find_nearest_neighbors(
    query: str,
    targets: list[str],
    n: int = 5,
    fp_type: str = "morgan",
) -> list[tuple[str, float]]:
    """Find the n nearest neighbors to a query compound.

    Args:
        query: Query SMILES string.
        targets: List of candidate SMILES strings.
        n: Number of neighbors to return.
        fp_type: Fingerprint type.

    Returns:
        Top-n (smiles, similarity) tuples sorted by descending similarity.
    """
    ranked = bulk_tanimoto(query, targets, fp_type=fp_type)
    return ranked[:n]


# @MX:ANCHOR: Public API for SMARTS-based substructure searching
# @MX:REASON: Called by CompoundCollection.filter_by_substructure and external users
def substructure_search(smiles_list: list[str], smarts: str) -> list[bool]:
    """Check which SMILES match a SMARTS substructure pattern.

    Args:
        smiles_list: List of SMILES strings to search.
        smarts: SMARTS pattern to match.

    Returns:
        List of booleans, True if the compound matches the pattern.

    Raises:
        ValueError: If smarts is invalid or cannot be parsed.
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts!r}")

    results: list[bool] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append(False)
        else:
            results.append(mol.HasSubstructMatch(pattern))
    return results


def substructure_match_atoms(smiles: str, smarts: str) -> list[tuple[int, ...]]:
    """Get all matching atom index tuples for a SMARTS pattern in a molecule.

    Args:
        smiles: SMILES string of the molecule to search.
        smarts: SMARTS pattern to match.

    Returns:
        List of tuples, each tuple contains the atom indices for one match.
        Returns an empty list if the molecule does not match.

    Raises:
        ValueError: If smarts is invalid or smiles is invalid.
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts!r}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    return list(mol.GetSubstructMatches(pattern))
