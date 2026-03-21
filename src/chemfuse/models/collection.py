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

        Returns:
            New CompoundCollection with filtered compounds.
        """
        filtered: list[Compound] = []

        for compound in self.compounds:
            props = compound.properties

            if mw_range is not None and props.molecular_weight is not None:
                if not (mw_range[0] <= props.molecular_weight <= mw_range[1]):
                    continue

            if logp_range is not None and props.xlogp is not None:
                if not (logp_range[0] <= props.xlogp <= logp_range[1]):
                    continue

            if tpsa_range is not None and props.tpsa is not None:
                if not (tpsa_range[0] <= props.tpsa <= tpsa_range[1]):
                    continue

            if hbd_max is not None and props.hbd_count is not None:
                if props.hbd_count > hbd_max:
                    continue

            if hba_max is not None and props.hba_count is not None:
                if props.hba_count > hba_max:
                    continue

            if rotatable_bonds_max is not None and props.rotatable_bonds is not None:
                if props.rotatable_bonds > rotatable_bonds_max:
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
