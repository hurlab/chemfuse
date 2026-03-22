"""R-group decomposition and SAR table generation.

Provides R-group decomposition against a core scaffold (SMARTS) and
SAR table generation joining R-group data with bioactivity values.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

from chemfuse.core.exceptions import OptionalDependencyError

if TYPE_CHECKING:
    from chemfuse.models.compound import Compound

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional RDKit import
# ---------------------------------------------------------------------------

_Chem = None
_RGroupDecomposition = None
_RDKIT_AVAILABLE = False

try:
    from rdkit import Chem as _Chem  # type: ignore[import-not-found]
    from rdkit.Chem import rdRGroupDecomposition as _rdRGD  # type: ignore[import-not-found]

    _RGroupDecomposition = _rdRGD.RGroupDecomposition
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


def decompose_rgroups(smiles_list: list[str], core_smarts: str) -> pd.DataFrame:
    """Decompose a list of SMILES into R-groups against a core scaffold.

    Compounds that do not match the core scaffold are silently excluded.
    Invalid SMILES strings are silently excluded.

    Args:
        smiles_list: List of SMILES strings to decompose.
        core_smarts: SMARTS string defining the core scaffold.

    Returns:
        DataFrame with columns: ``smiles``, ``Core``, ``R1``, ``R2``, ...
        (R-group columns as SMILES strings). Empty DataFrame if no compounds
        match the core or if ``smiles_list`` is empty.

    Raises:
        ValueError: If ``core_smarts`` is not a valid SMARTS pattern.
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles_list:
        return pd.DataFrame()

    # Validate core SMARTS
    core = _Chem.MolFromSmarts(core_smarts)
    if core is None:
        raise ValueError(f"Invalid core SMARTS: {core_smarts!r}")

    # Parse SMILES, track original SMILES for each valid mol
    mols: list[object] = []
    original_smiles: list[str] = []
    for smi in smiles_list:
        mol = _Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            original_smiles.append(smi)

    if not mols:
        return pd.DataFrame()

    # Run R-group decomposition
    decomp = _RGroupDecomposition(core)

    matched_indices: list[int] = []
    for i, mol in enumerate(mols):
        # Add() returns the row index >= 0 on success, -1 on failure
        if decomp.Add(mol) >= 0:
            matched_indices.append(i)

    if not matched_indices:
        return pd.DataFrame()

    decomp.Process()

    # GetRGroupsAsColumns returns dict: {label: [Mol | None, ...]}
    groups: dict[str, list[object | None]] = decomp.GetRGroupsAsColumns()

    # Convert Mol objects to SMILES strings (None becomes empty string)
    groups_smiles: dict[str, list[str]] = {}
    for label, mol_list in groups.items():
        groups_smiles[label] = [
            _Chem.MolToSmiles(m) if m is not None else "" for m in mol_list
        ]

    # Build result DataFrame
    matched_smiles = [original_smiles[i] for i in matched_indices]
    result = {"smiles": matched_smiles}
    result.update(groups_smiles)

    return pd.DataFrame(result)


def rgroup_sar_table(
    compounds: list[Compound],
    core_smarts: str,
    activity_type: str = "IC50",
) -> pd.DataFrame:
    """Generate a SAR table with R-groups and bioactivity values.

    Decomposes compounds into R-groups against ``core_smarts``, then joins
    with bioactivity data filtered by ``activity_type``. Compounds without
    matching activity are included with NaN activity values.

    Args:
        compounds: List of Compound objects to process.
        core_smarts: SMARTS string defining the core scaffold.
        activity_type: Activity type to extract (e.g., ``"IC50"``, ``"Ki"``).
            Case-insensitive matching is applied.

    Returns:
        DataFrame with columns: ``smiles``, ``Core``, ``R1``, ``R2``, ...,
        ``activity_value``, ``activity_units``, ``pic50``.
        Compounds not matching the core are excluded.
        Compounds without activity data have NaN in activity columns.

    Raises:
        ValueError: If ``core_smarts`` is not a valid SMARTS pattern.
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not compounds:
        return pd.DataFrame()

    smiles_list = [c.smiles for c in compounds]

    # Run R-group decomposition
    rgroup_df = decompose_rgroups(smiles_list, core_smarts)

    if rgroup_df.empty:
        return pd.DataFrame()

    # Build a lookup from SMILES -> best activity for the given activity_type
    activity_map: dict[str, dict[str, float | str | None]] = {}
    activity_type_lower = activity_type.lower()

    for compound in compounds:
        smi = compound.smiles
        if not smi or smi in activity_map:
            continue

        best_value: float | None = None
        best_units: str | None = None
        best_pic50: float | None = None

        for bio in getattr(compound, "bioactivities", []):
            bio_type = getattr(bio, "activity_type", None)
            if bio_type is None:
                continue
            if bio_type.lower() != activity_type_lower:
                continue

            # Use value_nm for comparison to pick the best (lowest) value
            value_nm = getattr(bio, "value_nm", None)
            raw_value = getattr(bio, "value", None)

            if raw_value is None:
                continue

            # Select the entry with the smallest value_nm (strongest activity)
            if best_value is None:
                best_value = raw_value
                best_units = getattr(bio, "units", None)
                best_pic50 = getattr(bio, "pic50", None)
            elif value_nm is not None:
                # Compare using value_nm when available
                current_nm_val = None
                for prev_bio in compound.bioactivities:
                    if getattr(prev_bio, "value", None) == best_value:
                        current_nm_val = getattr(prev_bio, "value_nm", None)
                        break
                if current_nm_val is None or value_nm < current_nm_val:
                    best_value = raw_value
                    best_units = getattr(bio, "units", None)
                    best_pic50 = getattr(bio, "pic50", None)

        activity_map[smi] = {
            "activity_value": best_value,
            "activity_units": best_units,
            "pic50": best_pic50,
        }

    # Join activity data onto the R-group DataFrame
    rgroup_df["activity_value"] = rgroup_df["smiles"].map(
        lambda s: activity_map.get(s, {}).get("activity_value")
    )
    rgroup_df["activity_units"] = rgroup_df["smiles"].map(
        lambda s: activity_map.get(s, {}).get("activity_units")
    )
    rgroup_df["pic50"] = rgroup_df["smiles"].map(
        lambda s: activity_map.get(s, {}).get("pic50")
    )

    return rgroup_df


__all__ = [
    "decompose_rgroups",
    "rgroup_sar_table",
]
