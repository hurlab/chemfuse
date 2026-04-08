"""Pandas DataFrame accessor for ChemFuse.

Registers the `.cf` accessor on pandas DataFrames, enabling cheminformatics
operations directly on DataFrames that contain a SMILES column.

Usage:
    import pandas as pd
    import chemfuse  # registers the .cf accessor

    df = pd.DataFrame({"smiles": ["CC(=O)Oc1ccccc1C(=O)O"]})
    desc_df = df.cf.compute_descriptors()
    filtered = df.cf.filter_druglike("lipinski")
"""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional dependency guards
# ---------------------------------------------------------------------------

try:
    import pandas as pd

    _PANDAS_AVAILABLE = True
except ImportError:
    _PANDAS_AVAILABLE = False

try:
    from rdkit import Chem  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


def _require_rdkit(method: str) -> None:
    """Raise OptionalDependencyError if RDKit is not installed."""
    if not _RDKIT_AVAILABLE:
        from chemfuse.core.exceptions import OptionalDependencyError

        raise OptionalDependencyError("rdkit", "compute")


# ---------------------------------------------------------------------------
# SMILES column detection helpers
# ---------------------------------------------------------------------------

_SMILES_COL_NAMES = frozenset({"smiles", "canonical_smiles", "canonicalsmiles", "smi"})


def _detect_smiles_col(df: pd.DataFrame) -> str | None:  # type: ignore[name-defined]
    """Return the first column whose lower-cased name matches a known SMILES alias."""
    for col in df.columns:
        if str(col).lower() in _SMILES_COL_NAMES:
            return col
    return None


# ---------------------------------------------------------------------------
# Supported drug-likeness rules
# ---------------------------------------------------------------------------

_SUPPORTED_RULES = ("lipinski", "veber", "ghose", "egan", "muegge")


def _apply_druglike_rule(rule: str, smiles: str) -> bool:
    """Return True if the SMILES passes the named drug-likeness rule."""
    from chemfuse.compute.druglikeness import (
        egan_filter,
        ghose_filter,
        lipinski_filter,
        muegge_filter,
        veber_filter,
    )

    filter_fn = {
        "lipinski": lipinski_filter,
        "veber": veber_filter,
        "ghose": ghose_filter,
        "egan": egan_filter,
        "muegge": muegge_filter,
    }[rule]

    result = filter_fn({}, smiles=smiles)
    return result.pass_filter


# ---------------------------------------------------------------------------
# Accessor
# ---------------------------------------------------------------------------

if _PANDAS_AVAILABLE:

    @pd.api.extensions.register_dataframe_accessor("cf")
    class ChemFuseAccessor:
        """Pandas accessor for ChemFuse cheminformatics operations.

        Registered as ``df.cf``. All mutating methods return new DataFrames;
        the original is never modified.

        The accessor auto-detects the SMILES column by looking for column names
        ``smiles``, ``canonical_smiles``, ``canonicalsmiles``, or ``smi``
        (case-insensitive).

        Raises:
            ValueError: At construction time when no SMILES column is found.
        """

        def __init__(self, pandas_obj: pd.DataFrame) -> None:
            self._df = pandas_obj
            self._smiles_col = _detect_smiles_col(pandas_obj)

        # ------------------------------------------------------------------
        # Internal helpers
        # ------------------------------------------------------------------

        def _get_smiles_col(self) -> str:
            """Return the SMILES column name, raising ValueError if absent."""
            if self._smiles_col is None:
                raise ValueError(
                    "No SMILES column found. Expected one of: "
                    + ", ".join(sorted(_SMILES_COL_NAMES))
                    + f". DataFrame columns: {list(self._df.columns)}"
                )
            return self._smiles_col

        def _iter_smiles(self) -> list[tuple[int, Any, str | None]]:
            """Yield (original_index, row_label, smiles_or_None) tuples.

            Returns None for the SMILES value when the cell is NaN/None/empty
            or when the SMILES string is chemically invalid.
            """
            col = self._get_smiles_col()
            result = []
            for idx in self._df.index:
                raw = self._df.at[idx, col]
                if pd.isna(raw) or not str(raw).strip():
                    result.append((idx, idx, None))
                    continue
                smi = str(raw).strip()
                if _RDKIT_AVAILABLE:
                    mol = Chem.MolFromSmiles(smi)
                    if mol is None:
                        logger.warning("Invalid SMILES at index %s: %r", idx, smi)
                        result.append((idx, idx, None))
                        continue
                result.append((idx, idx, smi))
            return result

        # ------------------------------------------------------------------
        # Public API
        # ------------------------------------------------------------------

        def compute_descriptors(self) -> pd.DataFrame:
            """Compute RDKit molecular descriptors for each row.

            Adds 200+ descriptor columns to the DataFrame. Rows with invalid
            SMILES receive NaN for all descriptor columns.

            Returns:
                New DataFrame with original columns plus descriptor columns.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("compute_descriptors")
            from chemfuse.compute.descriptors import compute_descriptors

            records: list[dict[str, Any]] = []
            for idx, _, smi in self._iter_smiles():
                if smi is None:
                    records.append({})
                else:
                    try:
                        records.append(compute_descriptors(smi))
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("compute_descriptors failed at index %s: %s", idx, exc)
                        records.append({})

            desc_df = pd.DataFrame(records, index=self._df.index)
            return pd.concat([self._df, desc_df], axis=1)

        def compute_fingerprints(
            self,
            fp_type: str = "morgan",
            n_bits: int = 2048,
        ) -> pd.DataFrame:
            """Add a 'fingerprint' column containing bit-vector dictionaries.

            Each cell in the 'fingerprint' column contains the dict returned by
            ``compute_fingerprint()`` (keys: type, bits, num_bits, num_on_bits),
            or None for invalid SMILES.

            Args:
                fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit',
                    'topological_torsion', 'atom_pair').
                n_bits: Bit vector length (ignored for MACCS which is 167-bit).

            Returns:
                New DataFrame with an added 'fingerprint' column.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("compute_fingerprints")
            from chemfuse.compute.fingerprints import compute_fingerprint

            fps: list[Any] = []
            for idx, _, smi in self._iter_smiles():
                if smi is None:
                    fps.append(None)
                else:
                    try:
                        fps.append(compute_fingerprint(smi, fp_type=fp_type, n_bits=n_bits))
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("compute_fingerprint failed at index %s: %s", idx, exc)
                        fps.append(None)

            result = self._df.copy()
            result["fingerprint"] = fps
            return result

        def filter_druglike(self, rule: str = "lipinski") -> pd.DataFrame:
            """Filter rows to only drug-like compounds.

            Args:
                rule: Drug-likeness rule to apply. One of 'lipinski', 'veber',
                    'ghose', 'egan', 'muegge'.

            Returns:
                Filtered DataFrame containing only passing rows.

            Raises:
                ValueError: If no SMILES column is detected or the rule is
                    unknown.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("filter_druglike")
            rule_lower = rule.lower()
            if rule_lower not in _SUPPORTED_RULES:
                raise ValueError(
                    f"Unknown drug-likeness rule: {rule!r}. "
                    f"Supported: {', '.join(_SUPPORTED_RULES)}"
                )

            mask: list[bool] = []
            for _, _, smi in self._iter_smiles():
                if smi is None:
                    mask.append(False)
                else:
                    try:
                        mask.append(_apply_druglike_rule(rule_lower, smi))
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("filter_druglike failed: %s", exc)
                        mask.append(False)

            return self._df.loc[mask].copy()

        def predict_admet(self) -> pd.DataFrame:
            """Add ADMET prediction columns to the DataFrame.

            Adds one column per ADMET property (solubility, bbb_permeability,
            cyp1a2_inhibition, etc.) plus an 'admet_overall_score' column.
            Values are the numeric prediction value from each ADMETPrediction.
            Rows with invalid SMILES receive NaN.

            Returns:
                New DataFrame with original columns plus ADMET columns.

            Raises:
                ValueError: If no SMILES column is detected.
            """
            from chemfuse.compute.admet import predict_admet as _predict_admet

            records: list[dict[str, Any]] = []
            for idx, _, smi in self._iter_smiles():
                if smi is None:
                    records.append({})
                else:
                    try:
                        profile = _predict_admet(smi)
                        row: dict[str, Any] = {
                            f"admet_{name}": pred.value
                            for name, pred in profile.predictions.items()
                        }
                        row["admet_overall_score"] = profile.overall_score
                        records.append(row)
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("predict_admet failed at index %s: %s", idx, exc)
                        records.append({})

            admet_df = pd.DataFrame(records, index=self._df.index)
            return pd.concat([self._df, admet_df], axis=1)

        def add_scaffolds(self) -> pd.DataFrame:
            """Add 'scaffold' and 'generic_scaffold' columns.

            Uses Bemis-Murcko decomposition. Rows with invalid SMILES or
            acyclic molecules (no scaffold) receive None.

            Returns:
                New DataFrame with added 'scaffold' and 'generic_scaffold'
                columns.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("add_scaffolds")
            from chemfuse.analyze.scaffolds import generic_scaffold, murcko_scaffold

            scaffolds: list[str | None] = []
            generic_scaffolds: list[str | None] = []

            for _, _, smi in self._iter_smiles():
                if smi is None:
                    scaffolds.append(None)
                    generic_scaffolds.append(None)
                else:
                    scaffolds.append(murcko_scaffold(smi))
                    generic_scaffolds.append(generic_scaffold(smi))

            result = self._df.copy()
            result["scaffold"] = scaffolds
            result["generic_scaffold"] = generic_scaffolds
            return result

        def standardize(self) -> pd.DataFrame:
            """Add a 'standardized_smiles' column.

            Applies salt stripping and tautomer canonicalization. Rows with
            invalid SMILES receive None.

            Returns:
                New DataFrame with an added 'standardized_smiles' column.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("standardize")
            from chemfuse.compute.standardization import standardize_mol

            std_smiles: list[str | None] = []
            for _, _, smi in self._iter_smiles():
                if smi is None:
                    std_smiles.append(None)
                else:
                    try:
                        std_smiles.append(standardize_mol(smi))
                    except Exception as exc:  # noqa: BLE001
                        logger.warning("standardize failed: %s", exc)
                        std_smiles.append(None)

            result = self._df.copy()
            result["standardized_smiles"] = std_smiles
            return result

        def to_collection(self) -> Any:
            """Convert the DataFrame to a CompoundCollection.

            Maps columns to Compound fields using these rules:
            - The SMILES column maps to ``Compound.smiles``.
            - 'name' column maps to ``Compound.name``.
            - 'inchikey' / 'inchi_key' maps to ``Compound.inchikey``.
            - 'inchi' maps to ``Compound.inchi``.
            - 'cid' maps to ``Compound.cid`` (cast to int).
            - 'chembl_id' maps to ``Compound.chembl_id``.
            - 'formula' maps to ``Compound.formula``.

            Returns:
                CompoundCollection built from the DataFrame rows.

            Raises:
                ValueError: If no SMILES column is detected.
            """
            from chemfuse.models.collection import CompoundCollection
            from chemfuse.models.compound import Compound

            smiles_col = self._get_smiles_col()
            compounds: list[Compound] = []

            col_lower_map: dict[str, str] = {str(c).lower(): str(c) for c in self._df.columns}

            def _get(row: pd.Series, *aliases: str) -> Any:  # type: ignore[type-arg]
                for alias in aliases:
                    col = col_lower_map.get(alias)
                    if col is not None:
                        val = row.get(col)
                        if val is not None and not (
                            isinstance(val, float) and pd.isna(val)
                        ):
                            return val
                return None

            for idx in self._df.index:
                row = self._df.loc[idx]
                raw_smi = row.get(smiles_col, "")
                smi = "" if pd.isna(raw_smi) else str(raw_smi).strip()

                cid_val = _get(row, "cid")
                cid_int: int | None = None
                if cid_val is not None:
                    try:
                        cid_int = int(cid_val)
                    except (ValueError, TypeError):
                        pass

                compounds.append(
                    Compound(
                        smiles=smi,
                        name=_get(row, "name"),
                        inchikey=_get(row, "inchikey", "inchi_key"),
                        inchi=_get(row, "inchi"),
                        cid=cid_int,
                        chembl_id=_get(row, "chembl_id"),
                        formula=_get(row, "formula"),
                    )
                )

            return CompoundCollection(compounds=compounds)

        def diversity_pick(
            self,
            n: int,
            fp_type: str = "morgan",
        ) -> pd.DataFrame:
            """Select the n most diverse rows using MaxMin picking.

            Args:
                n: Number of diverse compounds to select.
                fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

            Returns:
                Subset DataFrame with the n most diverse rows.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("diversity_pick")
            from chemfuse.analyze.diversity import maxmin_pick

            smiles_list = [
                smi if smi is not None else ""
                for _, _, smi in self._iter_smiles()
            ]

            picks = maxmin_pick(smiles_list, n_pick=n, fp_type=fp_type)
            selected_labels = [self._df.index[i] for i in picks]
            return self._df.loc[selected_labels].copy()

        def tanimoto_matrix(
            self,
            fp_type: str = "morgan",
        ) -> Any:
            """Compute the pairwise Tanimoto similarity matrix.

            Rows/columns with invalid SMILES produce NaN rows and columns.

            Args:
                fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

            Returns:
                numpy ndarray of shape (n, n) with Tanimoto similarities.

            Raises:
                ValueError: If no SMILES column is detected.
                OptionalDependencyError: If RDKit is not installed.
            """
            _require_rdkit("tanimoto_matrix")
            from chemfuse.analyze.similarity import tanimoto_matrix as _tanimoto_matrix

            smiles_list = [
                smi if smi is not None else ""
                for _, _, smi in self._iter_smiles()
            ]
            return _tanimoto_matrix(smiles_list, fp_type=fp_type)

else:
    # pandas is not installed; define a placeholder so imports do not crash.
    class ChemFuseAccessor:  # type: ignore[no-redef]
        """Placeholder accessor when pandas is not installed."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise ImportError(
                "pandas is required to use the ChemFuse DataFrame accessor. "
                "Install it with: pip install pandas"
            )
