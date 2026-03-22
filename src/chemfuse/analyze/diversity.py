"""Diversity picking and diversity metrics for compound collections."""

from __future__ import annotations

import logging

import numpy as np

from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem  # type: ignore[import-not-found]
    from rdkit.Chem import AllChem, MACCSkeys  # type: ignore[import-not-found]
    from rdkit.SimDivFilters import rdSimDivPickers  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


def _require_rdkit() -> None:
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "compute")


def _get_bitvec_fp(smiles: str, fp_type: str = "morgan") -> object:
    """Compute an RDKit bit-vector fingerprint for a SMILES string.

    Args:
        smiles: SMILES string.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

    Returns:
        RDKit ExplicitBitVect fingerprint object.

    Raises:
        ValueError: If the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    if fp_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    if fp_type == "rdkit":
        return Chem.RDKFingerprint(mol)
    # Default: Morgan radius=2, 2048 bits
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


# @MX:ANCHOR: Public API entry point for MaxMin diversity picking
# @MX:REASON: Called by CompoundCollection.pick_diverse and external library design workflows
def maxmin_pick(
    smiles_list: list[str],
    n_pick: int,
    fp_type: str = "morgan",
) -> list[int]:
    """Select n maximally diverse compounds using the MaxMin algorithm.

    Uses RDKit's MaxMinPicker to select compounds that are maximally
    spread across chemical space. Invalid SMILES are silently excluded
    and the returned indices map back to valid entries in smiles_list.

    Args:
        smiles_list: List of SMILES strings.
        n_pick: Number of compounds to pick.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

    Returns:
        List of indices into smiles_list identifying the selected compounds.
        If n_pick >= len(smiles_list), all valid-SMILES indices are returned.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles_list:
        return []

    # Build fingerprints, tracking which indices produced valid fps
    valid_indices: list[int] = []
    fps: list[object] = []
    for idx, smi in enumerate(smiles_list):
        try:
            fp = _get_bitvec_fp(smi, fp_type)
            fps.append(fp)
            valid_indices.append(idx)
        except (ValueError, Exception) as exc:
            logger.warning("maxmin_pick: skipping index %d (%r): %s", idx, smi, exc)

    if not fps:
        return []

    n_valid = len(fps)
    if n_pick >= n_valid:
        return list(valid_indices)

    picker = rdSimDivPickers.MaxMinPicker()
    # LazyBitVectorPick expects a flat callable that returns a fingerprint by index
    picks = picker.LazyBitVectorPick(fps, n_valid, n_pick)

    # Map local indices back to original smiles_list indices
    return [valid_indices[i] for i in picks]


# @MX:ANCHOR: Public API entry point for collection diversity scoring
# @MX:REASON: Called by CompoundCollection.diversity_score and analysis workflows
def diversity_score(
    smiles_list: list[str],
    fp_type: str = "morgan",
) -> float:
    """Compute the mean pairwise Tanimoto distance of a compound collection.

    Diversity is defined as ``1 - mean_pairwise_tanimoto_similarity``.
    A score of 0.0 means all compounds are identical; ~1.0 means maximally
    diverse.

    Args:
        smiles_list: List of SMILES strings.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

    Returns:
        Mean pairwise Tanimoto distance in [0.0, 1.0].
        Returns 0.0 for collections with fewer than 2 compounds.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if len(smiles_list) < 2:
        return 0.0

    # Import here to avoid circular dependency at module load time
    from chemfuse.analyze.similarity import tanimoto_matrix

    matrix = tanimoto_matrix(smiles_list, fp_type=fp_type)

    n = matrix.shape[0]
    if n < 2:
        return 0.0

    # Extract upper-triangle (excluding diagonal) for pairwise similarities
    upper_indices = np.triu_indices(n, k=1)
    pairwise_sims = matrix[upper_indices]

    # Remove NaN values (from invalid SMILES)
    valid_sims = pairwise_sims[~np.isnan(pairwise_sims)]
    if len(valid_sims) == 0:
        return 0.0

    mean_similarity = float(np.mean(valid_sims))
    return 1.0 - mean_similarity
