"""ChEMBL REST API source adapter."""

from __future__ import annotations

import logging
from typing import Any
from urllib.parse import quote

from chemfuse.core.exceptions import NotFoundError
from chemfuse.models.bioactivity import Bioactivity
from chemfuse.models.compound import Compound, CompoundProperties
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_RATE_LIMIT = 1.0  # requests per second

logger = logging.getLogger(__name__)


def _safe_float(value: Any) -> float | None:
    """Convert a value to float, returning None on failure."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class ChEMBLAdapter(SourceAdapter):
    """ChEMBL REST API adapter.

    Implements search by compound name, SMILES string, or ChEMBL ID.
    Supports retrieval of molecule details, bioactivity data, and mechanism
    of action for approved drugs.

    Rate limit: ~1 req/sec (no API key required).
    API base: https://www.ebi.ac.uk/chembl/api/data/
    """

    name = "chembl"
    base_url = CHEMBL_BASE_URL
    rate_limit = CHEMBL_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the ChEMBL adapter.

        Args:
            cache: Optional Cache instance.
            timeout: HTTP request timeout in seconds.
        """
        self._http = AsyncHTTPClient(
            source_name=self.name,
            base_url=self.base_url,
            rate_limit=self.rate_limit,
            timeout=timeout,
            cache=cache,
        )

    async def search(
        self,
        query: str,
        query_type: str = "name",
    ) -> list[Compound]:
        """Search ChEMBL for compounds.

        Args:
            query: Search query string.
            query_type: One of 'name', 'smiles', 'identifier' (ChEMBL ID).

        Returns:
            List of Compound objects.

        Raises:
            NotFoundError: If no compounds match the query.
        """
        qt = query_type.lower()

        if qt in ("identifier", "chembl_id", "id") or (
            isinstance(query, str) and query.upper().startswith("CHEMBL")
        ):
            compound = await self.get_by_id(query)
            if compound is None:
                return []
            return [compound]

        if qt == "smiles":
            return await self._search_by_smiles(query)

        # Default: search by name
        return await self._search_by_name(query)

    async def _search_by_name(self, name: str) -> list[Compound]:
        """Search molecules by preferred name."""
        data = await self._http.get(
            "/molecule/search",
            params={"q": name, "format": "json"},
        )
        return self._parse_molecule_list(data)

    async def _search_by_smiles(self, smiles: str) -> list[Compound]:
        """Search molecules by SMILES string."""
        data = await self._http.get(
            "/molecule/search",
            params={"q": smiles, "format": "json"},
        )
        return self._parse_molecule_list(data)

    async def get_by_id(self, chembl_id: str) -> Compound | None:
        """Get a molecule by its ChEMBL ID.

        Args:
            chembl_id: ChEMBL ID (e.g., "CHEMBL25").

        Returns:
            Compound object, or None if not found.

        Raises:
            NotFoundError: If the compound is not found.
        """
        try:
            data = await self._http.get(
                f"/molecule/{quote(chembl_id, safe='')}",
                params={"format": "json"},
            )
        except NotFoundError:
            return None

        return self._parse_molecule(data)

    async def get_bioactivities(
        self,
        chembl_id: str,
        limit: int = 1000,
    ) -> list[Bioactivity]:
        """Retrieve bioactivity data for a compound.

        Handles pagination automatically: follows next URLs until all results
        are collected.

        Args:
            chembl_id: ChEMBL molecule ID.
            limit: Results per page (max 1000).

        Returns:
            List of Bioactivity objects.
        """
        activities: list[Bioactivity] = []
        url = "/activity"
        params: dict[str, Any] = {
            "molecule_chembl_id": quote(chembl_id, safe=""),
            "format": "json",
            "limit": min(limit, 1000),
        }

        max_pages = 100
        page_count = 0
        while True:
            page_count += 1
            if page_count > max_pages:
                logger.warning(
                    "ChEMBL pagination limit reached (%d pages)", max_pages,
                )
                break

            data = await self._http.get(url, params=params)
            if not isinstance(data, dict):
                break

            page_activities = data.get("activities", [])
            for act in page_activities:
                bioact = self._parse_activity(act, chembl_id)
                if bioact is not None:
                    activities.append(bioact)

            # Follow pagination
            page_meta = data.get("page_meta", {})
            next_url = page_meta.get("next") if page_meta else None
            if next_url and next_url.startswith("https://www.ebi.ac.uk/"):
                # next_url is a full URL; use it directly
                url = next_url
                params = None  # URL already includes query string
            else:
                break

        return activities

    async def get_mechanism(self, chembl_id: str) -> list[dict]:
        """Retrieve mechanism of action data for approved drugs.

        Args:
            chembl_id: ChEMBL molecule ID.

        Returns:
            List of mechanism dicts with target_name, action_type, references.
        """
        data = await self._http.get(
            "/mechanism",
            params={"molecule_chembl_id": chembl_id, "format": "json"},
        )
        if not isinstance(data, dict):
            return []

        mechanisms = data.get("mechanisms", [])
        result = []
        for mech in mechanisms:
            result.append({
                "target_name": mech.get("target_chembl_id", ""),
                "action_type": mech.get("action_type", ""),
                "mechanism_of_action": mech.get("mechanism_of_action", ""),
                "references": mech.get("mechanism_refs", []),
            })
        return result

    async def get_properties(self, identifier: str) -> dict:
        """Get detailed properties for a compound by ChEMBL ID.

        Args:
            identifier: ChEMBL ID.

        Returns:
            Dictionary of property names to values.
        """
        compound = await self.get_by_id(identifier)
        if compound is None:
            return {}
        return compound.to_dict()

    def is_available(self) -> bool:
        """Always returns True (availability checked at request time)."""
        return True

    # --- Private helpers ---

    def _parse_molecule_list(self, data: Any) -> list[Compound]:
        """Parse a molecules list response into Compound objects."""
        if not isinstance(data, dict):
            return []

        molecules = data.get("molecules", [])
        compounds = []
        for mol in molecules:
            compound = self._parse_molecule(mol)
            if compound is not None:
                compounds.append(compound)
        return compounds

    def _parse_molecule(self, mol: Any) -> Compound | None:
        """Parse a single molecule dict into a Compound.

        Args:
            mol: Raw molecule dict from ChEMBL API.

        Returns:
            Compound object, or None if mol data is unusable.
        """
        if not isinstance(mol, dict):
            return None

        chembl_id = mol.get("molecule_chembl_id", "")
        if not chembl_id:
            return None

        # Extract structures
        structures = mol.get("molecule_structures") or {}
        smiles = structures.get("canonical_smiles", "") or ""
        inchikey = structures.get("standard_inchi_key") or structures.get("inchi_key")

        # Extract properties
        props_raw = mol.get("molecule_properties") or {}
        properties = CompoundProperties(
            molecular_weight=_safe_float(props_raw.get("full_mwt")),
            xlogp=_safe_float(props_raw.get("alogp")),
            hbd_count=self._safe_int(props_raw.get("hbd")),
            hba_count=self._safe_int(props_raw.get("hba")),
            tpsa=_safe_float(props_raw.get("psa")),
            rotatable_bonds=self._safe_int(props_raw.get("rtb")),
            heavy_atom_count=self._safe_int(props_raw.get("heavy_atoms")),
        )

        name = mol.get("pref_name") or None

        # CF-E08: parse clinical/regulatory metadata from ChEMBL molecule data
        max_phase_raw = mol.get("max_phase")
        max_phase: int | None = None
        if max_phase_raw is not None:
            try:
                max_phase = int(max_phase_raw)
            except (TypeError, ValueError):
                pass

        molecule_type: str | None = mol.get("molecule_type") or None

        return Compound(
            chembl_id=chembl_id,
            smiles=smiles,
            inchikey=inchikey,
            name=name,
            formula=props_raw.get("molecular_formula") or mol.get("molecular_formula"),
            sources=["chembl"],
            properties=properties,
            max_phase=max_phase,
            molecule_type=molecule_type,
        )

    def _parse_activity(self, act: dict, molecule_chembl_id: str) -> Bioactivity | None:
        """Parse a single activity dict into a Bioactivity object.

        Args:
            act: Raw activity dict from ChEMBL API.
            molecule_chembl_id: The queried molecule's ChEMBL ID.

        Returns:
            Bioactivity object, or None if the activity data is incomplete.
        """
        target_name = act.get("target_pref_name") or act.get("target_chembl_id") or ""
        target_id = act.get("target_chembl_id")
        activity_type = act.get("standard_type") or act.get("type") or "other"
        units = act.get("standard_units") or act.get("units")
        relation = act.get("standard_relation") or act.get("relation")
        assay_type = act.get("assay_type")
        reference = act.get("document_chembl_id")

        value_raw = act.get("standard_value") or act.get("value")
        value = _safe_float(value_raw)

        assay_chembl_id = act.get("assay_chembl_id")

        # CF-E08: parse assay quality metadata
        confidence_score_raw = act.get("confidence_score")
        confidence_score: int | None = None
        if confidence_score_raw is not None:
            try:
                confidence_score = int(confidence_score_raw)
            except (TypeError, ValueError):
                pass
        assay_description: str | None = act.get("assay_description") or None
        data_validity_comment: str | None = act.get("data_validity_comment") or None

        return Bioactivity(
            target_name=target_name,
            target_id=target_id,
            activity_type=activity_type,
            value=value,
            units=units,
            relation=relation,
            assay_type=assay_type,
            source="chembl",
            reference=reference,
            confidence_score=confidence_score,
            assay_chembl_id=assay_chembl_id,
            assay_description=assay_description,
            data_validity_comment=data_validity_comment,
        )

    @staticmethod
    def _safe_int(value: Any) -> int | None:
        """Convert value to int, returning None on failure."""
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None


# backward compatibility alias
ChEMBLSource = ChEMBLAdapter

__all__ = ["ChEMBLAdapter", "ChEMBLSource"]
