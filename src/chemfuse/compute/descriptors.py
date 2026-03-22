"""RDKit-based molecular descriptor computation for ChemFuse.

Provides 200+ RDKit descriptors with graceful fallback when RDKit is not
installed. Invalid SMILES strings produce warnings and empty dictionaries
rather than exceptions during batch computation.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem  # type: ignore[import-not-found]
    from rdkit.Chem import Descriptors  # type: ignore[import-not-found]

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def compute_descriptors(smiles: str) -> dict[str, float]:
    """Compute 200+ molecular descriptors for a SMILES string using RDKit.

    Args:
        smiles: SMILES representation of the molecule.

    Returns:
        Dictionary mapping descriptor names (str) to values (float).
        Returns an empty dict if the SMILES cannot be parsed by RDKit.

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES string, returning empty descriptors: %r", smiles)
        return {}

    result: dict[str, float] = {}
    # Descriptors.descList is the canonical list of 200+ RDKit descriptors
    for name, fn in Descriptors.descList:
        try:
            value = fn(mol)
            if value is not None:
                result[name] = float(value)
        except Exception:  # noqa: BLE001
            # Some descriptors may fail for edge-case molecules; skip them
            pass

    return result


def compute_descriptors_batch(smiles_list: list[str]) -> list[dict[str, float]]:
    """Compute descriptors for multiple SMILES strings.

    Invalid SMILES produce an empty dict with a warning rather than raising.

    Args:
        smiles_list: List of SMILES strings.

    Returns:
        List of descriptor dictionaries (one per input SMILES).

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    results: list[dict[str, float]] = []
    for smi in smiles_list:
        results.append(compute_descriptors(smi))
    return results


def is_valid_smiles(smiles: str) -> bool:
    """Check whether a SMILES string is chemically valid.

    When RDKit is not available the function returns False because validity
    cannot be determined without a parser.

    Args:
        smiles: SMILES string to validate.

    Returns:
        True if the SMILES is valid, False if invalid or RDKit is not available.
    """
    if not RDKIT_AVAILABLE:
        return False

    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


def smiles_to_inchi(smiles: str) -> str | None:
    """Convert a SMILES string to an InChI string.

    Args:
        smiles: SMILES string.

    Returns:
        InChI string, or None if conversion fails.

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    from rdkit.Chem.inchi import MolToInchi  # type: ignore[import-not-found]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    result: str | None = MolToInchi(mol)
    return result


def smiles_to_inchikey(smiles: str) -> str | None:
    """Convert a SMILES string to an InChIKey.

    Args:
        smiles: SMILES string.

    Returns:
        InChIKey string, or None if conversion fails.

    Raises:
        ImportError: If RDKit is not installed.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("Install rdkit: pip install chemfuse[rdkit]")

    from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi  # type: ignore[import-not-found]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchi = MolToInchi(mol)
    if inchi is None:
        return None
    result: str | None = InchiToInchiKey(inchi)
    return result
