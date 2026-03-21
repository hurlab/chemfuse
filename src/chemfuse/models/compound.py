"""Unified Compound data model for ChemFuse."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

if TYPE_CHECKING:
    from chemfuse.models.prediction import DrugLikeness

logger = logging.getLogger(__name__)


class CompoundProperties(BaseModel):
    """Physicochemical properties of a compound."""

    model_config = ConfigDict(populate_by_name=True)

    molecular_weight: float | None = None
    exact_mass: float | None = None
    xlogp: float | None = None
    tpsa: float | None = None
    hbd_count: int | None = None  # hydrogen bond donors
    hba_count: int | None = None  # hydrogen bond acceptors
    rotatable_bonds: int | None = None
    heavy_atom_count: int | None = None
    aromatic_rings: int | None = None
    complexity: float | None = None


class Compound(BaseModel):
    """Unified compound representation from any chemical database source.

    Tracks identifiers, properties, and which sources contributed data.
    """

    model_config = ConfigDict(populate_by_name=True)

    # Primary identifiers
    cid: int | None = None                    # PubChem CID
    chembl_id: str | None = None              # ChEMBL ID
    smiles: str = ""                           # Canonical SMILES (required)
    inchi: str | None = None                   # InChI string
    inchikey: str | None = None                # InChIKey
    name: str | None = None                    # Common name
    formula: str | None = None                 # Molecular formula
    synonyms: list[str] = Field(default_factory=list)

    # Source tracking
    sources: list[str] = Field(default_factory=list)

    # Physicochemical properties
    properties: CompoundProperties = Field(default_factory=CompoundProperties)

    # Bioactivity and binding data (populated via enrich())
    bioactivities: list[Any] = Field(default_factory=list)  # list[Bioactivity]
    binding_data: list[Any] = Field(default_factory=list)   # list[BindingMeasurement]
    target_associations: list[Any] = Field(default_factory=list)  # list[TargetAssociation]
    patents: list[Any] = Field(default_factory=list)              # list[Patent]

    # Computed fields (populated by compute methods)
    descriptors: dict[str, float] = Field(default_factory=dict)
    fingerprints: dict[str, Any] = Field(default_factory=dict)
    druglikeness: Any | None = None  # DrugLikeness | None

    def to_dict(self) -> dict[str, Any]:
        """Convert compound to a flat dictionary.

        Returns:
            Dictionary with all compound fields (including nested properties).
        """
        d = self.model_dump()
        # Flatten properties
        props = d.pop("properties", {})
        if props:
            for k, v in props.items():
                if v is not None:
                    d[k] = v
        return d

    def to_series(self) -> pd.Series:
        """Convert compound to a pandas Series (for use in DataFrames).

        Returns:
            pandas Series with compound data.
        """
        return pd.Series(self.to_dict())

    def merge(self, other: Compound) -> Compound:
        """Merge data from another Compound into this one.

        The other compound must represent the same molecule (matching InChIKey
        or CID). Non-None fields from other are merged into this compound.

        Args:
            other: Another Compound object to merge data from.

        Returns:
            New merged Compound.

        Raises:
            ValueError: If compounds do not represent the same molecule.
        """
        # Check if they represent the same molecule
        same = False
        if self.inchikey and other.inchikey and self.inchikey == other.inchikey:
            same = True
        elif self.cid and other.cid and self.cid == other.cid:
            same = True
        elif self.smiles and other.smiles and self.smiles == other.smiles:
            same = True

        if not same:
            if self.inchikey and other.inchikey and self.inchikey != other.inchikey:
                raise ValueError(
                    f"Cannot merge compounds with different InChIKeys: "
                    f"{self.inchikey!r} vs {other.inchikey!r}"
                )
            if self.smiles and other.smiles:
                # If we have SMILES but no other match, consider different molecules
                raise ValueError(
                    f"Cannot merge compounds with different identifiers: "
                    f"{self.inchikey!r} vs {other.inchikey!r}"
                )

        # Merge fields
        merged_data = self.model_dump()
        other_data = other.model_dump()

        for field in ["cid", "chembl_id", "inchi", "inchikey", "name", "formula"]:
            if merged_data.get(field) is None and other_data.get(field) is not None:
                merged_data[field] = other_data[field]

        # Merge SMILES (prefer non-empty)
        if not merged_data.get("smiles") and other_data.get("smiles"):
            merged_data["smiles"] = other_data["smiles"]

        # Merge synonyms
        existing_synonyms = set(merged_data.get("synonyms", []))
        for syn in other_data.get("synonyms", []):
            existing_synonyms.add(syn)
        merged_data["synonyms"] = list(existing_synonyms)

        # Merge sources
        merged_sources = list(set(merged_data.get("sources", [])) | set(other_data.get("sources", [])))
        merged_data["sources"] = merged_sources

        # Merge properties (other fills in missing values)
        self_props = merged_data.get("properties", {})
        other_props = other_data.get("properties", {})
        for k, v in other_props.items():
            if self_props.get(k) is None and v is not None:
                self_props[k] = v
        merged_data["properties"] = self_props

        # Merge bioactivities and binding_data (combine lists, deduplicate is left to caller)
        self_bioactivities = merged_data.get("bioactivities", [])
        other_bioactivities = other_data.get("bioactivities", [])
        merged_data["bioactivities"] = self_bioactivities + [
            b for b in other_bioactivities if b not in self_bioactivities
        ]

        self_binding = merged_data.get("binding_data", [])
        other_binding = other_data.get("binding_data", [])
        merged_data["binding_data"] = self_binding + [
            b for b in other_binding if b not in self_binding
        ]

        return Compound(**merged_data)

    async def enrich(
        self,
        sources: list[str] | None = None,
        binding: bool = False,
        targets: bool = False,
        patents: bool = False,
        force: bool = False,
    ) -> None:
        """Fetch additional data from specified sources and merge into this compound.

        Lazy loading: only explicitly requested data is fetched. Already-fetched
        sources are skipped unless force=True. Independent source fetches run in
        parallel via asyncio.gather with per-source error isolation.

        Args:
            sources: List of source names to fetch data from (e.g., ["chembl"]).
            binding: Whether to fetch binding data from BindingDB.
            targets: Whether to fetch disease-target associations from Open Targets.
            patents: Whether to fetch patent data from SureChEMBL.
            force: Re-fetch even if source already present in self.sources.
        """
        import asyncio

        from chemfuse.sources import registry

        if sources is None:
            sources = []

        for source_name in sources:
            if not force and source_name in self.sources:
                continue
            try:
                adapter = registry.get(source_name)
                # Fetch bioactivities from ChEMBL
                if source_name == "chembl":
                    get_bioactivities = getattr(adapter, "get_bioactivities", None)
                    if get_bioactivities is not None:
                        chembl_id = self.chembl_id
                        if not chembl_id:
                            # Try to look up by SMILES or name
                            search_query = self.inchikey or self.smiles or self.name
                            if search_query:
                                qt = "inchi" if self.inchikey else "smiles" if self.smiles else "name"
                                results = await adapter.search(search_query, query_type=qt)
                                if results and results[0].chembl_id:
                                    chembl_id = results[0].chembl_id
                                    if not self.chembl_id:
                                        self.chembl_id = chembl_id

                        if chembl_id:
                            bioacts = await get_bioactivities(chembl_id)
                            self.bioactivities = self.bioactivities + bioacts
                            if source_name not in self.sources:
                                self.sources = self.sources + [source_name]
            except Exception as exc:
                logger.warning("Enrichment from %s failed: %s", source_name, exc)

        # Build independent parallel tasks for binding, targets, patents
        tasks = []
        task_names = []

        if binding and (force or "bindingdb" not in self.sources):
            tasks.append(self._enrich_binding(registry))
            task_names.append("bindingdb")

        if targets and (force or "opentargets" not in self.sources):
            tasks.append(self._enrich_targets(registry))
            task_names.append("opentargets")

        if patents and (force or "surechembl" not in self.sources):
            tasks.append(self._enrich_patents(registry))
            task_names.append("surechembl")

        if tasks:
            results = await asyncio.gather(*tasks, return_exceptions=True)
            for task_name, result in zip(task_names, results, strict=True):
                if isinstance(result, Exception):
                    logger.warning("Enrichment from %s failed: %s", task_name, result)

    async def _enrich_binding(self, registry: Any) -> None:
        """Fetch binding data from BindingDB and merge into this compound."""
        adapter = registry.get("bindingdb")
        search_by_smiles = getattr(adapter, "search_by_smiles", None)
        if search_by_smiles is not None and self.smiles:
            measurements = await search_by_smiles(self.smiles)
            self.binding_data = self.binding_data + measurements
            if "bindingdb" not in self.sources:
                self.sources = self.sources + ["bindingdb"]

    async def _enrich_targets(self, registry: Any) -> None:
        """Fetch disease-target associations from Open Targets and merge into this compound.

        Resolves ChEMBL ID via UniChem if not directly available (R-ENGINE-04).
        """
        chembl_id = self.chembl_id

        # R-ENGINE-04: resolve ChEMBL ID if not present
        if not chembl_id:
            chembl_id = await self._resolve_chembl_id(registry)
            if chembl_id:
                self.chembl_id = chembl_id

        if not chembl_id:
            logger.warning(
                "Cannot enrich targets: no ChEMBL ID available for compound %s",
                self.name or self.smiles or "unknown",
            )
            return

        adapter = registry.get("opentargets")
        search_by_chembl_id = getattr(adapter, "search_by_chembl_id", None)
        if search_by_chembl_id is not None:
            associations = await search_by_chembl_id(chembl_id)
            self.target_associations = self.target_associations + associations
            if "opentargets" not in self.sources:
                self.sources = self.sources + ["opentargets"]

    async def _enrich_patents(self, registry: Any) -> None:
        """Fetch patent data from SureChEMBL and merge into this compound."""
        adapter = registry.get("surechembl")

        # Try InChIKey first, then SMILES
        if self.inchikey:
            search_by_inchikey = getattr(adapter, "search_by_inchikey", None)
            if search_by_inchikey is not None:
                patent_list = await search_by_inchikey(self.inchikey)
                self.patents = self.patents + patent_list
                if "surechembl" not in self.sources:
                    self.sources = self.sources + ["surechembl"]
                return

        if self.smiles:
            search_by_smiles = getattr(adapter, "search_by_smiles", None)
            if search_by_smiles is not None:
                patent_list = await search_by_smiles(self.smiles)
                self.patents = self.patents + patent_list
                if "surechembl" not in self.sources:
                    self.sources = self.sources + ["surechembl"]

    async def _resolve_chembl_id(self, registry: Any) -> str | None:
        """Attempt to resolve a ChEMBL ID via UniChem cross-reference.

        Resolution chain (SPEC Technical Notes):
        1. Check self.chembl_id (already done by caller)
        2. If InChIKey available, query UniChem by InChIKey
        3. If CID available, query UniChem by PubChem CID -> ChEMBL mapping
        4. If resolution fails, return None

        Args:
            registry: SourceRegistry instance.

        Returns:
            ChEMBL ID string or None.
        """
        try:
            unichem = registry.get("unichem")
            cross_reference = getattr(unichem, "cross_reference", None)
            map_identifiers = getattr(unichem, "map_identifiers", None)

            if self.inchikey and cross_reference is not None:
                mappings = await cross_reference(self.inchikey)
                chembl_id = mappings.get("chembl")
                if chembl_id:
                    return str(chembl_id)

            if self.cid and map_identifiers is not None:
                mappings = await map_identifiers(str(self.cid), "pubchem")
                chembl_id = mappings.get("chembl")
                if chembl_id:
                    return str(chembl_id)

        except Exception as exc:
            logger.debug("ChEMBL ID resolution via UniChem failed: %s", exc)

        return None

    # ------------------------------------------------------------------
    # Computation methods (SPEC-CF-002)
    # ------------------------------------------------------------------

    def compute_descriptors(self) -> dict[str, float]:
        """Compute RDKit molecular descriptors and store in self.descriptors.

        Returns:
            Dictionary of descriptor name to value.

        Raises:
            ImportError: If RDKit is not installed.
            ValueError: If this compound has no SMILES.
        """
        from chemfuse.compute.descriptors import compute_descriptors as _compute

        if not self.smiles:
            raise ValueError("Compound has no SMILES string; cannot compute descriptors.")

        self.descriptors = _compute(self.smiles)
        return self.descriptors

    def compute_fingerprints(
        self,
        types: list[str] | None = None,
        radius: int = 2,
        n_bits: int = 2048,
    ) -> dict[str, Any]:
        """Compute molecular fingerprints and store in self.fingerprints.

        Args:
            types: List of fingerprint type names to compute.  Defaults to
                ['morgan', 'maccs'].
            radius: Radius for Morgan fingerprints.
            n_bits: Bit vector length (not used for MACCS).

        Returns:
            Dictionary mapping fingerprint type to result dict.

        Raises:
            ImportError: If RDKit is not installed.
        """
        from chemfuse.compute.fingerprints import compute_fingerprint as _compute_fp

        if types is None:
            types = ["morgan", "maccs"]

        if not self.smiles:
            raise ValueError("Compound has no SMILES string; cannot compute fingerprints.")

        for fp_type in types:
            result = _compute_fp(self.smiles, fp_type=fp_type, radius=radius, n_bits=n_bits)
            if result is not None:
                self.fingerprints[fp_type] = result

        return self.fingerprints

    def check_drug_likeness(self) -> DrugLikeness:
        """Evaluate drug-likeness using all five standard filters.

        Uses pre-populated properties when available; falls back to RDKit
        computation from SMILES when available.

        Returns:
            DrugLikeness with individual FilterResult for each filter.
        """
        from chemfuse.compute.druglikeness import check_drug_likeness as _check

        # Build a properties dict from self.properties (PubChem-style names)
        props = self.properties.model_dump()
        # Normalise names to what druglikeness.py understands
        props_norm: dict[str, Any] = {
            "molecular_weight": props.get("molecular_weight"),
            "xlogp": props.get("xlogp"),
            "tpsa": props.get("tpsa"),
            "hbd_count": props.get("hbd_count"),
            "hba_count": props.get("hba_count"),
            "rotatable_bonds": props.get("rotatable_bonds"),
            "heavy_atom_count": props.get("heavy_atom_count"),
        }

        result = _check(props_norm, smiles=self.smiles or None)
        self.druglikeness = result
        return result

    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        name = self.name or self.smiles[:20] if self.smiles else "unknown"
        cid_str = f", cid={self.cid}" if self.cid else ""
        return f"Compound(name={name!r}{cid_str})"
