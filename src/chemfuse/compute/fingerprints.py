"""Molecular fingerprint generation for ChemFuse.

Supports Morgan (ECFP), MACCS, RDKit, TopologicalTorsion, and AtomPair
fingerprints. Gracefully degrades when RDKit is not installed.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem, DataStructs  # type: ignore[import-not-found]
    from rdkit.Chem import (  # type: ignore[import-not-found]
        AllChem,
        MACCSkeys,
        rdMolDescriptors,  # type: ignore[import-not-found]
    )

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Supported fingerprint type identifiers
_FP_MORGAN = "morgan"
_FP_MACCS = "maccs"
_FP_RDKIT = "rdkit"
_FP_TOPOLOGICAL_TORSION = "topological_torsion"
_FP_ATOM_PAIR = "atom_pair"

SUPPORTED_FP_TYPES = (
    _FP_MORGAN,
    _FP_MACCS,
    _FP_RDKIT,
    _FP_TOPOLOGICAL_TORSION,
    _FP_ATOM_PAIR,
)


def fp_matrix(smiles_list: list[str], fp_type: str = "morgan", n_bits: int = 2048) -> np.ndarray:
    """Build a numpy bit-vector matrix from a list of SMILES strings.

    Handles MACCS (167-bit), RDKit, and Morgan fingerprints.  Invalid SMILES
    produce an all-zero row rather than raising an exception.

    Args:
        smiles_list: List of SMILES strings.
        fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').
        n_bits: Bit vector length for morgan/rdkit fingerprints.

    Returns:
        numpy ndarray of shape (len(smiles_list), actual_bits).

    Raises:
        ImportError: If RDKit is not installed.
    """
    import logging

    import numpy as np

    _logger = logging.getLogger(__name__)

    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    # Determine actual bit length from fp_type
    if fp_type == "maccs":
        actual_bits = 167
    elif fp_type == "rdkit":
        actual_bits = n_bits
    else:
        actual_bits = n_bits

    rows = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                rows.append(np.zeros(actual_bits, dtype=np.uint8))
                continue
            if fp_type == "maccs":
                fp = MACCSkeys.GenMACCSKeys(mol)
                arr = np.zeros(167, dtype=np.uint8)
            elif fp_type == "rdkit":
                fp = Chem.RDKFingerprint(mol)
                arr = np.zeros(n_bits, dtype=np.uint8)
            else:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
                arr = np.zeros(n_bits, dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            rows.append(arr)
        except Exception as exc:
            _logger.warning("fp_matrix: skipping %r: %s", smi, exc)
            rows.append(np.zeros(actual_bits, dtype=np.uint8))
    return np.array(rows, dtype=float)


def compute_fingerprint(
    smiles: str,
    fp_type: str = _FP_MORGAN,
    radius: int = 2,
    n_bits: int = 2048,
) -> dict[str, Any] | None:
    """Generate a molecular fingerprint for a SMILES string.

    Args:
        smiles: SMILES representation of the molecule.
        fp_type: Fingerprint type. One of 'morgan', 'maccs', 'rdkit',
            'topological_torsion', 'atom_pair'.
        radius: Radius for Morgan fingerprints (default 2 = ECFP4).
        n_bits: Bit vector length (ignored for MACCS which is always 167-bit).

    Returns:
        Dictionary with keys:
            - 'type': fingerprint type name
            - 'bits': list of set bit positions
            - 'num_bits': total number of bits
            - 'num_on_bits': number of set bits
        Returns None if the SMILES cannot be parsed (logs a warning).

    Raises:
        ImportError: If RDKit is not installed.
        ValueError: If fp_type is not a supported type.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    fp_type_lower = fp_type.lower()
    if fp_type_lower not in SUPPORTED_FP_TYPES:
        raise ValueError(
            f"Unknown fingerprint type: {fp_type!r}. "
            f"Supported: {', '.join(SUPPORTED_FP_TYPES)}"
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES string, returning None fingerprint: %r", smiles)
        return None

    if fp_type_lower == _FP_MORGAN:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    elif fp_type_lower == _FP_MACCS:
        fp = MACCSkeys.GenMACCSKeys(mol)
        n_bits = fp.GetNumBits()
    elif fp_type_lower == _FP_RDKIT:
        fp = Chem.RDKFingerprint(mol, fpSize=n_bits)
    elif fp_type_lower == _FP_TOPOLOGICAL_TORSION:
        fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol, nBits=n_bits
        )
    else:  # atom_pair
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol, nBits=n_bits
        )

    on_bits = list(fp.GetOnBits())
    return {
        "type": fp_type_lower,
        "bits": on_bits,
        "num_bits": fp.GetNumBits(),
        "num_on_bits": len(on_bits),
    }


def tanimoto_similarity(
    smiles_a: str,
    smiles_b: str,
    fp_type: str = _FP_MORGAN,
    radius: int = 2,
    n_bits: int = 2048,
) -> float:
    """Compute Tanimoto similarity between two molecules.

    Args:
        smiles_a: SMILES of the first molecule.
        smiles_b: SMILES of the second molecule.
        fp_type: Fingerprint type to use for similarity.
        radius: Radius for Morgan fingerprints.
        n_bits: Bit vector length.

    Returns:
        Tanimoto similarity coefficient (0.0 to 1.0).

    Raises:
        ImportError: If RDKit is not installed.
        ValueError: If either SMILES string is invalid.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)

    if mol_a is None:
        raise ValueError(f"Invalid SMILES string: {smiles_a!r}")
    if mol_b is None:
        raise ValueError(f"Invalid SMILES string: {smiles_b!r}")

    fp_type_lower = fp_type.lower()

    if fp_type_lower == _FP_MORGAN:
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=n_bits)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=n_bits)
    elif fp_type_lower == _FP_MACCS:
        fp_a = MACCSkeys.GenMACCSKeys(mol_a)
        fp_b = MACCSkeys.GenMACCSKeys(mol_b)
    elif fp_type_lower == _FP_RDKIT:
        fp_a = Chem.RDKFingerprint(mol_a, fpSize=n_bits)
        fp_b = Chem.RDKFingerprint(mol_b, fpSize=n_bits)
    elif fp_type_lower == _FP_TOPOLOGICAL_TORSION:
        fp_a = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol_a, nBits=n_bits
        )
        fp_b = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol_b, nBits=n_bits
        )
    elif fp_type_lower == _FP_ATOM_PAIR:
        fp_a = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol_a, nBits=n_bits
        )
        fp_b = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol_b, nBits=n_bits
        )
    else:
        # Fallback to Morgan for unrecognised type
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=n_bits)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=n_bits)

    similarity: float = DataStructs.TanimotoSimilarity(fp_a, fp_b)
    return similarity


def bulk_tanimoto(
    query_smiles: str,
    target_smiles_list: list[str],
    fp_type: str = _FP_MORGAN,
    radius: int = 2,
    n_bits: int = 2048,
) -> list[tuple[str, float]]:
    """Compute Tanimoto similarity between a query and multiple target molecules.

    Invalid target SMILES strings are silently skipped.

    Args:
        query_smiles: SMILES of the query molecule.
        target_smiles_list: List of target SMILES strings.
        fp_type: Fingerprint type to use.
        radius: Radius for Morgan fingerprints.
        n_bits: Bit vector length.

    Returns:
        List of (target_smiles, similarity) tuples, sorted descending by similarity.

    Raises:
        ImportError: If RDKit is not installed.
        ValueError: If the query SMILES is invalid.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid query SMILES: {query_smiles!r}")

    fp_type_lower = fp_type.lower()

    # Build query fingerprint
    if fp_type_lower == _FP_MACCS:
        query_fp = MACCSkeys.GenMACCSKeys(query_mol)
    elif fp_type_lower == _FP_RDKIT:
        query_fp = Chem.RDKFingerprint(query_mol, fpSize=n_bits)
    elif fp_type_lower == _FP_TOPOLOGICAL_TORSION:
        query_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            query_mol, nBits=n_bits
        )
    elif fp_type_lower == _FP_ATOM_PAIR:
        query_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            query_mol, nBits=n_bits
        )
    else:
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius, nBits=n_bits)

    results: list[tuple[str, float]] = []
    for target_smi in target_smiles_list:
        target_mol = Chem.MolFromSmiles(target_smi)
        if target_mol is None:
            continue
        # Build target fingerprint using same type as query
        if fp_type_lower == _FP_MACCS:
            target_fp = MACCSkeys.GenMACCSKeys(target_mol)
        elif fp_type_lower == _FP_RDKIT:
            target_fp = Chem.RDKFingerprint(target_mol, fpSize=n_bits)
        elif fp_type_lower == _FP_TOPOLOGICAL_TORSION:
            target_fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                target_mol, nBits=n_bits
            )
        elif fp_type_lower == _FP_ATOM_PAIR:
            target_fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                target_mol, nBits=n_bits
            )
        else:
            target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, radius, nBits=n_bits)

        sim: float = DataStructs.TanimotoSimilarity(query_fp, target_fp)
        results.append((target_smi, sim))

    results.sort(key=lambda x: x[1], reverse=True)
    return results
