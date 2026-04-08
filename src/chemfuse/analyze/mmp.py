"""Matched Molecular Pair (MMP) analysis for ChemFuse.

Identifies pairs of compounds that differ by a single structural transformation
and correlates those transformations with changes in activity or properties.
"""

from __future__ import annotations

import logging
from itertools import combinations
from typing import Any

import pandas as pd

from chemfuse.core.exceptions import OptionalDependencyError

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem  # noqa: I001
    from rdkit.Chem import rdFMCS
    from rdkit.Chem import rdFingerprintGenerator

    _RDKIT_AVAILABLE = True
    _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
except ImportError:
    _RDKIT_AVAILABLE = False


def _require_rdkit() -> None:
    """Raise OptionalDependencyError when RDKit is not installed."""
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", "compute")


def _mol_from_smiles(smiles: str) -> Any | None:
    """Return an RDKit Mol object or None on failure."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except Exception:  # noqa: BLE001
        return None


def _heavy_atom_count(mol: Any) -> int:
    """Return the number of heavy atoms in an RDKit Mol object."""
    return mol.GetNumHeavyAtoms()


def _tanimoto(mol_a: Any, mol_b: Any) -> float:
    """Return the Tanimoto similarity between two RDKit Mol objects."""
    from rdkit import DataStructs

    fp_a = _MORGAN_GEN.GetFingerprint(mol_a)
    fp_b = _MORGAN_GEN.GetFingerprint(mol_b)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


def _mcs_heavy_atoms(mol_a: Any, mol_b: Any) -> int:
    """Return the heavy-atom count of the Maximum Common Substructure.

    Returns 0 when the MCS search times out or finds no common structure.
    """
    result = rdFMCS.FindMCS(
        [mol_a, mol_b],
        timeout=5,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        completeRingsOnly=True,
        matchValences=True,
    )
    if result is None or result.canceled:
        return 0
    # result.numAtoms counts heavy atoms in the MCS SMARTS pattern
    return result.numAtoms


def _extract_transform(mol_a: Any, mol_b: Any, mcs_smarts: str) -> tuple[str, str]:
    """Extract the transformation fragment (non-MCS part) for each molecule.

    Returns (transform_a, transform_b) as SMILES strings.
    Returns ("*", "*") when extraction fails.
    """
    try:
        core_query = Chem.MolFromSmarts(mcs_smarts)
        if core_query is None:
            return ("*", "*")

        def _get_transform(mol: Any) -> str:
            match = mol.GetSubstructMatch(core_query)
            if not match:
                return "*"
            core_atom_indices = set(match)
            edit_mol = Chem.RWMol(mol)
            # Remove core atoms; what remains is the transform
            atoms_to_remove = sorted(core_atom_indices, reverse=True)
            for idx in atoms_to_remove:
                edit_mol.RemoveAtom(idx)
            frag = edit_mol.GetMol()
            if frag.GetNumHeavyAtoms() == 0:
                return "H"  # Hydrogen replacement (no substituent)
            return Chem.MolToSmiles(frag)

        return (_get_transform(mol_a), _get_transform(mol_b))
    except Exception:  # noqa: BLE001
        return ("*", "*")


def find_matched_pairs(
    smiles_list: list[str],
    activities: list[float] | None = None,
    min_heavy_atoms_core: int = 5,
    max_heavy_atoms_transform: int = 13,
) -> pd.DataFrame:
    """Find matched molecular pairs in a set of compounds.

    A matched molecular pair is two compounds that differ by a single
    structural transformation (one R-group change on a shared core).

    Args:
        smiles_list: SMILES strings for compounds.
        activities: Optional activity values (e.g., pIC50) for each compound.
        min_heavy_atoms_core: Minimum heavy atoms in the shared core.
        max_heavy_atoms_transform: Maximum heavy atoms in the transformation.

    Returns:
        DataFrame with columns: smiles_a, smiles_b, core_smarts, transform_a,
        transform_b, activity_a, activity_b, activity_delta (if activities
        provided). Sorted by |activity_delta| descending when activities are given,
        otherwise by (smiles_a, smiles_b).

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    _require_rdkit()

    if not smiles_list:
        return pd.DataFrame(
            columns=[
                "smiles_a", "smiles_b", "core_smarts",
                "transform_a", "transform_b",
                "activity_a", "activity_b", "activity_delta",
            ]
        )

    if len(smiles_list) < 2:
        return pd.DataFrame(
            columns=[
                "smiles_a", "smiles_b", "core_smarts",
                "transform_a", "transform_b",
                "activity_a", "activity_b", "activity_delta",
            ]
        )

    if activities is not None and len(activities) != len(smiles_list):
        raise ValueError(
            f"activities length ({len(activities)}) must match "
            f"smiles_list length ({len(smiles_list)})"
        )

    # Parse molecules; track valid indices
    mols: list[Any | None] = [_mol_from_smiles(s) for s in smiles_list]
    valid_indices = [i for i, m in enumerate(mols) if m is not None]

    if len(valid_indices) < 2:
        return pd.DataFrame(
            columns=[
                "smiles_a", "smiles_b", "core_smarts",
                "transform_a", "transform_b",
                "activity_a", "activity_b", "activity_delta",
            ]
        )

    rows: list[dict[str, Any]] = []

    # For large sets, pre-filter by Tanimoto similarity before MCS
    use_similarity_prefilter = len(valid_indices) > 100
    sim_prefilter_threshold = 0.5

    for idx_a, idx_b in combinations(valid_indices, 2):
        mol_a = mols[idx_a]
        mol_b = mols[idx_b]

        # Similarity pre-filter for performance with large compound sets
        if use_similarity_prefilter:
            sim = _tanimoto(mol_a, mol_b)
            if sim < sim_prefilter_threshold:
                continue

        ha_a = _heavy_atom_count(mol_a)
        ha_b = _heavy_atom_count(mol_b)

        mcs_n = _mcs_heavy_atoms(mol_a, mol_b)
        if mcs_n < min_heavy_atoms_core:
            continue

        # Transformation size check: each molecule must differ by at most
        # max_heavy_atoms_transform atoms from the MCS
        diff_a = ha_a - mcs_n
        diff_b = ha_b - mcs_n
        if diff_a > max_heavy_atoms_transform or diff_b > max_heavy_atoms_transform:
            continue

        # Retrieve the MCS SMARTS for transform extraction
        mcs_result = rdFMCS.FindMCS(
            [mol_a, mol_b],
            timeout=5,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            completeRingsOnly=True,
            matchValences=True,
        )
        if mcs_result is None or mcs_result.canceled or not mcs_result.smartsString:
            continue

        transform_a, transform_b = _extract_transform(
            mol_a, mol_b, mcs_result.smartsString
        )

        row: dict[str, Any] = {
            "smiles_a": smiles_list[idx_a],
            "smiles_b": smiles_list[idx_b],
            "core_smarts": mcs_result.smartsString,
            "transform_a": transform_a,
            "transform_b": transform_b,
            "activity_a": None,
            "activity_b": None,
            "activity_delta": None,
        }

        if activities is not None:
            act_a = activities[idx_a]
            act_b = activities[idx_b]
            row["activity_a"] = act_a
            row["activity_b"] = act_b
            row["activity_delta"] = act_b - act_a

        rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=[
                "smiles_a", "smiles_b", "core_smarts",
                "transform_a", "transform_b",
                "activity_a", "activity_b", "activity_delta",
            ]
        )

    df = pd.DataFrame(rows)

    if activities is not None and "activity_delta" in df.columns:
        df = df.sort_values(
            by="activity_delta",
            key=lambda col: col.abs(),
            ascending=False,
        ).reset_index(drop=True)

    return df


def summarize_transformations(
    pairs_df: pd.DataFrame,
    min_count: int = 2,
) -> pd.DataFrame:
    """Summarize recurring structural transformations and their average activity effect.

    Groups matched pairs by transformation type and computes statistics.

    Args:
        pairs_df: DataFrame from find_matched_pairs().
        min_count: Minimum number of occurrences for a transformation to be included.

    Returns:
        DataFrame with columns: transform_from, transform_to, count,
        mean_activity_delta, std_activity_delta, compounds_involved.
        Empty DataFrame when pairs_df is empty or has no recurring transformations.
    """
    if pairs_df.empty:
        return pd.DataFrame(
            columns=[
                "transform_from", "transform_to", "count",
                "mean_activity_delta", "std_activity_delta", "compounds_involved",
            ]
        )

    required_cols = {"transform_a", "transform_b"}
    if not required_cols.issubset(pairs_df.columns):
        return pd.DataFrame(
            columns=[
                "transform_from", "transform_to", "count",
                "mean_activity_delta", "std_activity_delta", "compounds_involved",
            ]
        )

    # Normalize transform direction: sort lexicographically so (A->B) and (B->A) group together
    def _canonical_transform(row: pd.Series) -> tuple[str, str]:
        ta, tb = str(row["transform_a"]), str(row["transform_b"])
        return (ta, tb) if ta <= tb else (tb, ta)

    pairs_df = pairs_df.copy()
    pairs_df[["transform_from", "transform_to"]] = pd.DataFrame(
        pairs_df.apply(_canonical_transform, axis=1).tolist(),
        index=pairs_df.index,
    )

    has_activity = (
        "activity_delta" in pairs_df.columns
        and pairs_df["activity_delta"].notna().any()
    )

    group_keys = ["transform_from", "transform_to"]
    grouped = pairs_df.groupby(group_keys)

    records: list[dict[str, Any]] = []
    for (tf, tt), group in grouped:
        count = len(group)
        if count < min_count:
            continue

        compounds = sorted(
            set(group["smiles_a"].tolist() + group["smiles_b"].tolist())
        )

        row: dict[str, Any] = {
            "transform_from": tf,
            "transform_to": tt,
            "count": count,
            "mean_activity_delta": None,
            "std_activity_delta": None,
            "compounds_involved": compounds,
        }

        if has_activity:
            deltas = group["activity_delta"].dropna()
            if not deltas.empty:
                row["mean_activity_delta"] = float(deltas.mean())
                row["std_activity_delta"] = float(deltas.std()) if len(deltas) > 1 else 0.0

        records.append(row)

    if not records:
        return pd.DataFrame(
            columns=[
                "transform_from", "transform_to", "count",
                "mean_activity_delta", "std_activity_delta", "compounds_involved",
            ]
        )

    result_df = pd.DataFrame(records).sort_values("count", ascending=False).reset_index(drop=True)
    return result_df
