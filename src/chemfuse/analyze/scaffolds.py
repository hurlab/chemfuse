"""Murcko scaffold decomposition and frequency analysis.

Provides Bemis-Murcko scaffold extraction, generic scaffold generation,
and scaffold frequency/grouping utilities for compound collections.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import TYPE_CHECKING

from chemfuse.core.exceptions import OptionalDependencyError

if TYPE_CHECKING:
    from chemfuse.models.compound import Compound

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional RDKit import
# ---------------------------------------------------------------------------

_Chem = None
_MurckoScaffold = None
_RDKIT_AVAILABLE = False

try:
    from rdkit import Chem as _Chem  # type: ignore[import-not-found]
    from rdkit.Chem.Scaffolds import (
        MurckoScaffold as _MurckoScaffold,  # type: ignore[import-not-found]
    )

    _RDKIT_AVAILABLE = True
except ImportError:
    pass


def _require_rdkit() -> None:
    """Raise OptionalDependencyError if RDKit is not installed."""
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "analyze")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def murcko_scaffold(smiles: str) -> str | None:
    """Extract the Bemis-Murcko scaffold from a SMILES string.

    Args:
        smiles: SMILES representation of the molecule.

    Returns:
        Canonical SMILES of the Murcko scaffold, or None if extraction fails
        (invalid SMILES, no scaffold, etc.).

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles:
        return None

    try:
        mol = _Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        scaffold_mol = _MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is None:
            return None

        scaffold_smiles = _Chem.MolToSmiles(scaffold_mol)
        # Empty string means the molecule has no ring/scaffold framework
        return scaffold_smiles if scaffold_smiles else None
    except Exception as exc:  # noqa: BLE001
        logger.debug("murcko_scaffold failed for %r: %s", smiles, exc)
        return None


def generic_scaffold(smiles: str) -> str | None:
    """Extract the generic (carbon-skeleton) Murcko scaffold from a SMILES string.

    The generic scaffold replaces all heteroatoms with carbon and removes
    bond-order information, producing a carbon-only ring framework.

    Args:
        smiles: SMILES representation of the molecule.

    Returns:
        Canonical SMILES of the generic scaffold, or None if extraction fails.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles:
        return None

    try:
        mol = _Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        scaffold_mol = _MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is None:
            return None

        generic_mol = _MurckoScaffold.MakeScaffoldGeneric(scaffold_mol)
        if generic_mol is None:
            return None

        generic_smiles = _Chem.MolToSmiles(generic_mol)
        return generic_smiles if generic_smiles else None
    except Exception as exc:  # noqa: BLE001
        logger.debug("generic_scaffold failed for %r: %s", smiles, exc)
        return None


def scaffold_frequency(smiles_list: list[str]) -> dict[str, int]:
    """Compute Murcko scaffold frequency across a list of SMILES.

    Scaffolds are extracted for each SMILES; invalid/None scaffolds are skipped.
    The result is sorted by frequency descending (most common scaffold first).

    Args:
        smiles_list: List of SMILES strings.

    Returns:
        Ordered dict mapping scaffold SMILES to occurrence count, sorted by
        count descending. Returns an empty dict for an empty input.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles_list:
        return {}

    counts: Counter[str] = Counter()
    for smi in smiles_list:
        scaffold = murcko_scaffold(smi)
        if scaffold is not None:
            counts[scaffold] += 1

    # Return sorted by frequency descending
    return dict(counts.most_common())


def group_by_scaffold(compounds: list[Compound]) -> dict[str, list[Compound]]:
    """Group Compound objects by their Murcko scaffold SMILES.

    Compounds whose scaffold cannot be extracted are placed under the key
    ``""`` (empty string).

    Args:
        compounds: List of Compound objects. Each must have a non-empty
            ``smiles`` attribute.

    Returns:
        Dict mapping scaffold SMILES (or ``""`` for unclassified) to a list
        of Compound objects sharing that scaffold.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    groups: dict[str, list[Compound]] = {}
    for compound in compounds:
        if compound.smiles:
            scaffold = murcko_scaffold(compound.smiles)
        else:
            scaffold = None

        key = scaffold if scaffold is not None else ""
        groups.setdefault(key, []).append(compound)

    return groups
