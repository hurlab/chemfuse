"""CompoundCollection model for managing sets of compounds."""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Callable, Iterator
from datetime import UTC, datetime
from typing import Any

import pandas as pd
from pydantic import BaseModel, Field

from chemfuse.models.compound import Compound

logger = logging.getLogger(__name__)


class CompoundCollection(BaseModel):
    """An ordered collection of Compound objects with metadata.

    Supports iteration, indexing, filtering, sorting, and export.
    """

    compounds: list[Compound] = Field(default_factory=list)
    query: str = ""
    sources: list[str] = Field(default_factory=list)
    timestamp: datetime = Field(default_factory=lambda: datetime.now(UTC))
    warnings: list[str] = Field(default_factory=list)

    def __len__(self) -> int:
        return len(self.compounds)

    def __iter__(self) -> Iterator[Compound]:
        return iter(self.compounds)

    def __getitem__(self, index: int | slice) -> Compound | list[Compound]:
        return self.compounds[index]

    def __repr__(self) -> str:
        return (
            f"CompoundCollection(n={len(self.compounds)}, "
            f"query={self.query!r}, sources={self.sources})"
        )

    def filter(
        self,
        mw_range: tuple[float, float] | None = None,
        logp_range: tuple[float, float] | None = None,
        tpsa_range: tuple[float, float] | None = None,
        hbd_max: int | None = None,
        hba_max: int | None = None,
        rotatable_bonds_max: int | None = None,
        sources: list[str] | None = None,
        strict: bool = False,
    ) -> CompoundCollection:
        """Filter compounds by property ranges.

        Args:
            mw_range: (min, max) molecular weight range.
            logp_range: (min, max) XLogP range.
            tpsa_range: (min, max) TPSA range.
            hbd_max: Maximum hydrogen bond donor count.
            hba_max: Maximum hydrogen bond acceptor count.
            rotatable_bonds_max: Maximum rotatable bond count.
            sources: Filter to compounds from specified sources only.
            strict: When True, compounds with None for any filtered property
                are excluded. When False (default), None values pass through.

        Returns:
            New CompoundCollection with filtered compounds.
        """
        filtered: list[Compound] = []

        for compound in self.compounds:
            props = compound.properties

            if mw_range is not None:
                if props.molecular_weight is None:
                    if strict:
                        continue
                elif not (mw_range[0] <= props.molecular_weight <= mw_range[1]):
                    continue

            if logp_range is not None:
                if props.xlogp is None:
                    if strict:
                        continue
                elif not (logp_range[0] <= props.xlogp <= logp_range[1]):
                    continue

            if tpsa_range is not None:
                if props.tpsa is None:
                    if strict:
                        continue
                elif not (tpsa_range[0] <= props.tpsa <= tpsa_range[1]):
                    continue

            if hbd_max is not None:
                if props.hbd_count is None:
                    if strict:
                        continue
                elif props.hbd_count > hbd_max:
                    continue

            if hba_max is not None:
                if props.hba_count is None:
                    if strict:
                        continue
                elif props.hba_count > hba_max:
                    continue

            if rotatable_bonds_max is not None:
                if props.rotatable_bonds is None:
                    if strict:
                        continue
                elif props.rotatable_bonds > rotatable_bonds_max:
                    continue

            if sources is not None:
                if not any(s in compound.sources for s in sources):
                    continue

            filtered.append(compound)

        return CompoundCollection(
            compounds=filtered,
            query=self.query,
            sources=self.sources,
            timestamp=self.timestamp,
        )

    def filter_druglike(
        self,
        rule: str = "lipinski",
    ) -> CompoundCollection:
        """Filter compounds by drug-likeness rules.

        Args:
            rule: Drug-likeness rule to apply. Options: 'lipinski', 'veber'.

        Returns:
            New CompoundCollection with drug-like compounds.
        """
        if rule == "lipinski":
            return self.filter(
                mw_range=(0, 500),
                logp_range=(-5, 5),
                hbd_max=5,
                hba_max=10,
            )
        elif rule == "veber":
            return self.filter(
                tpsa_range=(0, 140),
                rotatable_bonds_max=10,
            )
        else:
            return self

    def compute_all(
        self,
        descriptors: bool = True,
        fingerprints: bool = False,
        druglikeness: bool = True,
        progress_callback: Any | None = None,
    ) -> None:
        """Apply selected computations to all compounds in the collection.

        Computations that fail for a given compound (e.g. invalid SMILES) are
        skipped with a warning rather than aborting the batch.

        Args:
            descriptors: Compute RDKit molecular descriptors.
            fingerprints: Compute Morgan and MACCS fingerprints.
            druglikeness: Evaluate the five drug-likeness filters.
            progress_callback: Optional callable(completed: int, total: int,
                compound_name: str).
        """
        total = len(self.compounds)
        for idx, compound in enumerate(self.compounds):
            compound_name = compound.name or compound.smiles or str(compound.cid) or "unknown"
            try:
                if descriptors:
                    try:
                        compound.compute_descriptors()
                    except Exception as exc:  # noqa: BLE001
                        logger.warning(
                            "compute_all: descriptor computation failed for %s: %s",
                            compound_name, exc,
                        )
                if fingerprints:
                    try:
                        compound.compute_fingerprints()
                    except Exception as exc:  # noqa: BLE001
                        logger.warning(
                            "compute_all: fingerprint computation failed for %s: %s",
                            compound_name, exc,
                        )
                if druglikeness:
                    try:
                        compound.check_drug_likeness()
                    except Exception as exc:  # noqa: BLE001
                        logger.warning(
                            "compute_all: drug-likeness check failed for %s: %s",
                            compound_name, exc,
                        )
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "compute_all: unexpected error for %s: %s", compound_name, exc
                )

            if progress_callback is not None:
                progress_callback(idx + 1, total, compound_name)

    def filter_by_druglikeness(
        self,
        lipinski: bool = False,
        veber: bool = False,
        ghose: bool = False,
        egan: bool = False,
        muegge: bool = False,
    ) -> CompoundCollection:
        """Filter compounds by drug-likeness filter results.

        Compounds that do not yet have druglikeness computed will be evaluated
        on demand before filtering.

        Args:
            lipinski: Keep only compounds that pass Lipinski's Rule of Five.
            veber: Keep only compounds that pass the Veber filter.
            ghose: Keep only compounds that pass the Ghose filter.
            egan: Keep only compounds that pass the Egan filter.
            muegge: Keep only compounds that pass the Muegge filter.

        Returns:
            New CompoundCollection containing only compounds that pass all
            selected filters.
        """
        filtered: list[Compound] = []

        for compound in self.compounds:
            # Ensure drug-likeness is computed
            if compound.druglikeness is None:
                try:
                    compound.check_drug_likeness()
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "filter_by_druglikeness: cannot compute drug-likeness for %s: %s",
                        compound.name or compound.smiles, exc,
                    )
                    continue

            dl = compound.druglikeness
            if dl is None:
                continue

            if lipinski and not dl.lipinski.pass_filter:
                continue
            if veber and not dl.veber.pass_filter:
                continue
            if ghose and not dl.ghose.pass_filter:
                continue
            if egan and not dl.egan.pass_filter:
                continue
            if muegge and not dl.muegge.pass_filter:
                continue

            filtered.append(compound)

        return CompoundCollection(
            compounds=filtered,
            query=self.query,
            sources=self.sources,
            timestamp=self.timestamp,
        )

    def sort(
        self,
        by: str = "molecular_weight",
        ascending: bool = True,
    ) -> CompoundCollection:
        """Sort compounds by a property.

        Args:
            by: Property name to sort by. Supports compound fields and
                nested property fields (e.g., 'molecular_weight').
            ascending: Sort direction.

        Returns:
            New CompoundCollection with sorted compounds.
        """
        def get_sort_key(compound: Compound) -> Any:
            # Check compound-level fields first
            val = getattr(compound, by, None)
            if val is None:
                # Check properties
                val = getattr(compound.properties, by, None)
            # Put None values at the end
            return (val is None, val)

        sorted_compounds = sorted(
            self.compounds,
            key=get_sort_key,
            reverse=not ascending,
        )

        return CompoundCollection(
            compounds=sorted_compounds,
            query=self.query,
            sources=self.sources,
            timestamp=self.timestamp,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Convert collection to a pandas DataFrame.

        Returns:
            DataFrame with one row per compound and flattened properties.
        """
        if not self.compounds:
            return pd.DataFrame()

        records: list[dict[str, Any]] = []
        for compound in self.compounds:
            row = compound.to_dict()
            # Remove list fields for DataFrame simplicity
            row.pop("synonyms", None)
            row["sources"] = ", ".join(row.get("sources", []))
            records.append(row)

        df = pd.DataFrame(records)
        df = df.dropna(axis=1, how="all")
        return df

    async def enrich_all(
        self,
        sources: list[str] | None = None,
        binding: bool = False,
        force: bool = False,
        progress_callback: Callable[[int, int, str], None] | None = None,
    ) -> None:
        """Enrich all compounds with data from specified sources.

        Processes compounds sequentially respecting per-source rate limits.
        Failed individual enrichments are collected and do not abort the batch.

        Args:
            sources: Source names to fetch bioactivity data from.
            binding: Whether to fetch binding data from BindingDB.
            force: Re-fetch even if source already in compound.sources.
            progress_callback: Optional callable(completed, total, compound_name).
        """
        total = len(self.compounds)
        errors: list[tuple[str, Exception]] = []

        for idx, compound in enumerate(self.compounds):
            compound_name = compound.name or compound.smiles or str(compound.cid) or "unknown"
            try:
                await compound.enrich(sources=sources, binding=binding, force=force)
            except Exception as exc:
                errors.append((compound_name, exc))
                logger.warning("enrich_all: failed for %s: %s", compound_name, exc)

            if progress_callback is not None:
                progress_callback(idx + 1, total, compound_name)

            # Small delay between compounds to respect rate limits
            if idx < total - 1:
                await asyncio.sleep(0.1)

        if errors:
            logger.warning("enrich_all: %d compound(s) failed enrichment.", len(errors))

    def to_csv(self, output_path: str) -> str:
        """Export collection to CSV file.

        Args:
            output_path: Path for the output CSV file.

        Returns:
            Path to the created file.
        """
        from chemfuse.core.export import export_csv
        result = export_csv(self, output_path)
        return str(result)

    def to_excel(self, output_path: str, sheet_name: str = "ChemFuse Results") -> str:
        """Export collection to Excel file.

        Args:
            output_path: Path for the output Excel file.
            sheet_name: Excel sheet name.

        Returns:
            Path to the created file.
        """
        from chemfuse.core.export import export_excel
        result = export_excel(self, output_path, sheet_name=sheet_name)
        return str(result)

    def to_json(self, output_path: str) -> str:
        """Export collection to JSON file.

        Args:
            output_path: Path for the output JSON file.

        Returns:
            Path to the created file.
        """
        from chemfuse.core.export import export_json
        result = export_json(self, output_path)
        return str(result)

    # ------------------------------------------------------------------
    # Analysis methods (SPEC-CF-005)
    # ------------------------------------------------------------------

    def predict_admet(self) -> None:
        """Predict ADMET properties for all compounds in the collection.

        Results are stored in ``compound.admet_profile`` (a dict on the
        compound's extra data dict, keyed ``'admet'``).  The actual
        ADMETProfile is available as ``compound.admet``.
        """
        from chemfuse.compute.admet import predict_admet as _predict

        for compound in self.compounds:
            if not compound.smiles:
                continue
            try:
                profile = _predict(compound.smiles)
                # Store on compound as an attribute if possible, otherwise
                # use the extra dict pattern.
                object.__setattr__(compound, "_admet_profile", profile)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "predict_admet: failed for %s: %s",
                    compound.name or compound.smiles, exc,
                )

    def cluster(
        self,
        method: str = "butina",
        cutoff: float = 0.4,
        n_clusters: int = 5,
        fp_type: str = "morgan",
    ) -> list[int]:
        """Cluster compounds and attach cluster labels.

        Cluster labels are stored as ``compound.cluster_label`` (via
        ``_cluster_label`` attribute on each compound).

        Args:
            method: 'butina' or 'kmeans'.
            cutoff: Butina distance cutoff.
            n_clusters: Number of clusters for KMeans.
            fp_type: Fingerprint type.

        Returns:
            List of integer cluster labels.
        """
        from chemfuse.analyze.clustering import butina_clustering, kmeans_clustering

        smiles_list = [c.smiles for c in self.compounds]
        if method == "kmeans":
            labels = kmeans_clustering(smiles_list, n_clusters=n_clusters, fp_type=fp_type)
        else:
            labels = butina_clustering(smiles_list, cutoff=cutoff, fp_type=fp_type)

        for compound, label in zip(self.compounds, labels, strict=False):
            object.__setattr__(compound, "_cluster_label", label)

        return labels

    def reduce_dimensions(
        self,
        method: str = "umap",
        fp_type: str = "morgan",
    ) -> object:
        """Reduce compound fingerprints to 2D coordinates.

        Args:
            method: 'umap', 'tsne', or 'pca'.
            fp_type: Fingerprint type.

        Returns:
            numpy array of shape (n, 2).
        """

        from chemfuse.analyze.chemspace import reduce_dimensions as _reduce

        smiles_list = [c.smiles for c in self.compounds]
        return _reduce(smiles_list, method=method, fp_type=fp_type)

    def visualize_chemical_space(
        self,
        method: str = "umap",
        color_by: str = "cluster",
        fp_type: str = "morgan",
    ) -> object:
        """Visualize the chemical space of the collection.

        Args:
            method: Dimensionality reduction method.
            color_by: 'cluster' or a property name.
            fp_type: Fingerprint type.

        Returns:
            plotly.graph_objects.Figure.
        """

        from chemfuse.analyze.chemspace import plot_chemical_space
        from chemfuse.analyze.chemspace import reduce_dimensions as _reduce

        smiles_list = [c.smiles for c in self.compounds]
        coords = _reduce(smiles_list, method=method, fp_type=fp_type)

        labels = [c.name or c.smiles for c in self.compounds]
        colors: list[Any] | None = None

        if color_by == "cluster":
            cluster_labels = [
                getattr(c, "_cluster_label", None) for c in self.compounds
            ]
            if any(x is not None for x in cluster_labels):
                colors = [str(x) if x is not None else "?" for x in cluster_labels]

        return plot_chemical_space(coords, labels=labels, colors=colors)

    def filter_by_substructure(self, smarts: str) -> CompoundCollection:
        """Return a new collection containing only compounds matching the SMARTS pattern.

        Compounds with missing or invalid SMILES are excluded silently.

        Args:
            smarts: SMARTS substructure pattern to match.

        Returns:
            New CompoundCollection with only matching compounds.

        Raises:
            ValueError: If smarts is an invalid SMARTS pattern.
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.similarity import substructure_search

        smiles_list = [c.smiles or "" for c in self.compounds]
        # substructure_search returns False for invalid/missing SMILES
        matches = substructure_search(smiles_list, smarts)

        filtered = [
            compound
            for compound, matched in zip(self.compounds, matches, strict=False)
            if matched
        ]

        return CompoundCollection(
            compounds=filtered,
            query=self.query,
            sources=self.sources,
            timestamp=self.timestamp,
        )

    def decompose_rgroups(self, core_smarts: str) -> pd.DataFrame:
        """Decompose compounds into R-groups against a core scaffold.

        Compounds that do not match the core scaffold (SMARTS) are silently
        excluded. Requires RDKit to be installed.

        Args:
            core_smarts: SMARTS string defining the core scaffold.

        Returns:
            DataFrame with columns: ``smiles``, ``Core``, ``R1``, ``R2``, ...
            Empty DataFrame if no compounds match the core.

        Raises:
            ValueError: If ``core_smarts`` is not a valid SMARTS pattern.
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.rgroup import decompose_rgroups as _decompose

        smiles_list = [c.smiles for c in self.compounds]
        return _decompose(smiles_list, core_smarts)

    def sar_table(
        self,
        core_smarts: str,
        activity_type: str = "IC50",
    ) -> pd.DataFrame:
        """Generate a SAR table with R-groups and activity values.

        Decomposes compounds into R-groups against ``core_smarts``, then joins
        with bioactivity data filtered by ``activity_type``. Compounds without
        matching activity are included with NaN activity values.

        Args:
            core_smarts: SMARTS string defining the core scaffold.
            activity_type: Activity type to extract (e.g., ``"IC50"``, ``"Ki"``).

        Returns:
            DataFrame with columns: ``smiles``, ``Core``, ``R1``, ``R2``, ...,
            ``activity_value``, ``activity_units``, ``pic50``.

        Raises:
            ValueError: If ``core_smarts`` is not a valid SMARTS pattern.
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.rgroup import rgroup_sar_table as _sar_table

        return _sar_table(self.compounds, core_smarts, activity_type=activity_type)

    def detect_activity_cliffs(
        self,
        activity_col: str = "ic50",
        sim_threshold: float = 0.8,
        fp_type: str = "morgan",
    ) -> list[dict]:
        """Detect activity cliffs among compounds in the collection.

        Args:
            activity_col: Name of the activity column in compound.extra.
            sim_threshold: Minimum Tanimoto similarity.
            fp_type: Fingerprint type.

        Returns:
            List of cliff dicts sorted by descending cliff_score.
        """
        from chemfuse.analyze.sar import detect_activity_cliffs as _detect

        compound_dicts = []
        for c in self.compounds:
            d: dict[str, Any] = {"smiles": c.smiles}
            if c.name:
                d["name"] = c.name
            # Try to find the activity column in the compound's dict
            raw = c.to_dict()
            if activity_col in raw:
                d[activity_col] = raw[activity_col]
            elif hasattr(c, "bioactivities") and c.bioactivities:
                # Pull from first bioactivity record
                for bio in c.bioactivities:
                    if hasattr(bio, "value") and bio.value is not None:
                        d[activity_col] = bio.value
                        break
            compound_dicts.append(d)

        return _detect(compound_dicts, activity_col=activity_col, sim_threshold=sim_threshold, fp_type=fp_type)

    # ------------------------------------------------------------------
    # Scaffold methods (CF-E01)
    # ------------------------------------------------------------------

    def scaffold_frequency(self) -> pd.DataFrame:
        """Compute Murcko scaffold frequency across the collection.

        Scaffolds are computed on the fly from each compound's SMILES.
        Compounds without a valid SMILES or scaffold are excluded.

        Returns:
            DataFrame with columns ``scaffold``, ``count``, and ``percentage``,
            sorted by count descending. Returns an empty DataFrame when the
            collection is empty.

        Raises:
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.scaffolds import scaffold_frequency as _scaffold_freq

        smiles_list = [c.smiles for c in self.compounds if c.smiles]
        freq_dict = _scaffold_freq(smiles_list)

        if not freq_dict:
            return pd.DataFrame(columns=["scaffold", "count", "percentage"])

        total = sum(freq_dict.values())
        records = [
            {
                "scaffold": scaffold,
                "count": count,
                "percentage": round(count / total * 100, 2) if total > 0 else 0.0,
            }
            for scaffold, count in freq_dict.items()
        ]
        return pd.DataFrame(records)

    def group_by_scaffold(self) -> dict[str, CompoundCollection]:
        """Group compounds in this collection by their Murcko scaffold.

        Scaffolds are computed on the fly from each compound's SMILES.
        Compounds whose scaffold cannot be determined are grouped under the
        key ``""`` (empty string).

        Returns:
            Dict mapping scaffold SMILES to a :class:`CompoundCollection`
            containing all compounds that share that scaffold.

        Raises:
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.scaffolds import group_by_scaffold as _group

        groups = _group(self.compounds)

        return {
            scaffold_smi: CompoundCollection(
                compounds=compound_list,
                query=self.query,
                sources=self.sources,
                timestamp=self.timestamp,
            )
            for scaffold_smi, compound_list in groups.items()
        }

    # ------------------------------------------------------------------
    # Diversity methods (CF-E06)
    # ------------------------------------------------------------------

    def pick_diverse(
        self,
        n: int,
        fp_type: str = "morgan",
    ) -> CompoundCollection:
        """Select n maximally diverse compounds using the MaxMin algorithm.

        Args:
            n: Number of compounds to select.
            fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

        Returns:
            New CompoundCollection containing the n most diverse compounds.
            If n >= len(self), returns a copy of the full collection.

        Raises:
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.diversity import maxmin_pick as _maxmin_pick

        smiles_list = [c.smiles or "" for c in self.compounds]
        indices = _maxmin_pick(smiles_list, n_pick=n, fp_type=fp_type)
        picked = [self.compounds[i] for i in indices]

        return CompoundCollection(
            compounds=picked,
            query=self.query,
            sources=self.sources,
            timestamp=self.timestamp,
        )

    def diversity_score(self, fp_type: str = "morgan") -> float:
        """Compute mean pairwise Tanimoto distance (0=identical, ~1=diverse).

        Args:
            fp_type: Fingerprint type ('morgan', 'maccs', 'rdkit').

        Returns:
            Mean pairwise Tanimoto distance in [0.0, 1.0].
            Returns 0.0 for collections with fewer than 2 compounds.

        Raises:
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.diversity import diversity_score as _diversity_score

        smiles_list = [c.smiles for c in self.compounds if c.smiles]
        return _diversity_score(smiles_list, fp_type=fp_type)

    # ------------------------------------------------------------------
    # SAR Landscape methods (CF-E10)
    # ------------------------------------------------------------------

    def plot_sar_landscape(
        self,
        activity_col: str = "pic50",
        method: str = "umap",
        fp_type: str = "morgan",
        show_cliffs: bool = True,
        cliff_threshold: float = 0.8,
    ) -> dict[str, Any]:
        """Generate SAR landscape from compound bioactivity data.

        Extracts SMILES and activity values from the collection, then delegates
        to :func:`chemfuse.analyze.landscape.sar_landscape` to compute 2D
        coordinates, detect activity cliffs, and build an interactive figure.

        Activity values are sourced in the following priority order:

        1. ``pic50`` attribute on a :class:`~chemfuse.models.bioactivity.Bioactivity`
           record attached to the compound (when *activity_col* is ``"pic50"``).
        2. The first ``value_nm`` from a bioactivity record when pic50 is
           unavailable.
        3. A matching key in the compound's ``extra`` dict.

        Compounds without a resolvable activity value are silently skipped.

        Args:
            activity_col: Bioactivity attribute or extra-dict key to use.
                Defaults to ``"pic50"``.
            method: Dimensionality reduction method: "umap", "tsne", or "pca".
            fp_type: Fingerprint type for similarity and DR ("morgan", "maccs",
                "rdkit").
            show_cliffs: When *True*, cliff pairs are detected and included in
                the returned dict and drawn on the figure.
            cliff_threshold: Minimum Tanimoto similarity for cliff detection.

        Returns:
            dict as returned by :func:`~chemfuse.analyze.landscape.sar_landscape`.

        Raises:
            ValueError: If no compounds with valid activities are found.
        """
        from chemfuse.analyze.landscape import sar_landscape as _sar_landscape

        smiles_valid: list[str] = []
        activities_valid: list[float] = []

        for compound in self.compounds:
            if not compound.smiles:
                continue

            activity_value: float | None = None

            # Try bioactivities list
            if hasattr(compound, "bioactivities") and compound.bioactivities:
                for bio in compound.bioactivities:
                    if activity_col == "pic50" and hasattr(bio, "pic50") and bio.pic50 is not None:
                        activity_value = float(bio.pic50)
                        break
                    if activity_col == "value_nm" and hasattr(bio, "value_nm") and bio.value_nm is not None:
                        activity_value = float(bio.value_nm)
                        break
                    # Fallback: use first available value_nm as surrogate
                    if activity_value is None and hasattr(bio, "value_nm") and bio.value_nm is not None:
                        activity_value = float(bio.value_nm)

            # Fallback: try compound.extra dict
            if activity_value is None:
                raw = compound.to_dict()
                if activity_col in raw and raw[activity_col] is not None:
                    try:
                        activity_value = float(raw[activity_col])
                    except (TypeError, ValueError):
                        pass

            if activity_value is not None:
                smiles_valid.append(compound.smiles)
                activities_valid.append(activity_value)

        if not smiles_valid:
            raise ValueError(
                f"No compounds with valid '{activity_col}' activity values found in the collection. "
                "Ensure compounds have bioactivity data (SPEC-CF-002/CF-003)."
            )

        cliff_thresh = cliff_threshold if show_cliffs else 1.1  # > 1.0 means no cliffs detected

        return _sar_landscape(
            smiles_list=smiles_valid,
            activities=activities_valid,
            method=method,
            fp_type=fp_type,
            cliff_threshold=cliff_thresh,
        )

    # ------------------------------------------------------------------
    # MMP methods (Task 2)
    # ------------------------------------------------------------------

    def generate_sar_report(
        self,
        activity_type: str = "IC50",
        target: str | None = None,
        format: str = "markdown",
    ) -> str:
        """Generate a SAR report. See chemfuse.analyze.report for details.

        Args:
            activity_type: Activity type to analyze (IC50, Ki, etc.).
            target: Optional target name to filter bioactivities.
            format: Output format ("markdown" or "text").

        Returns:
            Formatted SAR report as a string.
        """
        from chemfuse.analyze.report import generate_sar_report as _gen
        return _gen(self, activity_type=activity_type, target=target, format=format)

    def find_matched_pairs(
        self,
        activity_type: str = "IC50",
        min_heavy_atoms_core: int = 5,
        max_heavy_atoms_transform: int = 13,
    ) -> pd.DataFrame:
        """Find matched molecular pairs in this collection.

        A matched molecular pair is two compounds that differ by a single
        structural transformation (one R-group change on a shared core).
        Activity values are extracted from the first bioactivity record that
        matches ``activity_type``.

        Args:
            activity_type: Activity type to extract (e.g., ``"IC50"``, ``"Ki"``).
                Used to look up activity values from compound bioactivity data.
            min_heavy_atoms_core: Minimum heavy atoms required in the shared core.
            max_heavy_atoms_transform: Maximum heavy atoms in the transformation.

        Returns:
            DataFrame with columns: smiles_a, smiles_b, core_smarts,
            transform_a, transform_b, activity_a, activity_b, activity_delta.
            activity_delta is None when no bioactivity data is available.

        Raises:
            OptionalDependencyError: If RDKit is not installed.
        """
        from chemfuse.analyze.mmp import find_matched_pairs as _find_mmp

        smiles_list = [c.smiles for c in self.compounds if c.smiles]
        compound_index = [i for i, c in enumerate(self.compounds) if c.smiles]

        activities: list[float] | None = None
        activity_values: list[float | None] = []

        for idx in compound_index:
            compound = self.compounds[idx]
            value: float | None = None
            if hasattr(compound, "bioactivities") and compound.bioactivities:
                for bio in compound.bioactivities:
                    bio_type = getattr(bio, "activity_type", None) or getattr(bio, "type", None)
                    if bio_type and activity_type.lower() in str(bio_type).lower():
                        v = getattr(bio, "value", None) or getattr(bio, "value_nm", None)
                        if v is not None:
                            try:
                                value = float(v)
                                break
                            except (TypeError, ValueError):
                                pass
            activity_values.append(value)

        # Only pass activities list when at least some values are available
        if any(v is not None for v in activity_values):
            # Replace None with NaN-safe float; MMP module handles None gracefully
            activities = [v if v is not None else float("nan") for v in activity_values]

        return _find_mmp(
            smiles_list=smiles_list,
            activities=activities,
            min_heavy_atoms_core=min_heavy_atoms_core,
            max_heavy_atoms_transform=max_heavy_atoms_transform,
        )
