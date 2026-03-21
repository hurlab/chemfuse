"""PubChem PUG-REST source adapter."""

from __future__ import annotations

import asyncio
import time
from typing import Any
from urllib.parse import quote

from chemfuse.core.exceptions import NotFoundError, SourceError, TimeoutError
from chemfuse.models.compound import Compound, CompoundProperties
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_RATE_LIMIT = 5.0  # requests per second

PROPERTY_LIST = ",".join([
    "MolecularFormula",
    "MolecularWeight",
    "ExactMass",
    "CanonicalSMILES",
    "IsomericSMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "HeavyAtomCount",
    "Complexity",
])


class PubChemAdapter(SourceAdapter):
    """PubChem PUG-REST API adapter.

    Implements search by name, SMILES, CID, formula, and InChI.
    Supports similarity search, substructure search, and bioassay retrieval.
    Rate limit: 5 requests/second.
    """

    name = "pubchem"
    base_url = PUBCHEM_BASE_URL
    rate_limit = PUBCHEM_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the PubChem adapter.

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
        """Search PubChem for compounds.

        Args:
            query: Search query string.
            query_type: One of 'name', 'smiles', 'cid', 'formula', 'inchi', 'identifier'.

        Returns:
            List of Compound objects.
        """
        qt = query_type.lower()

        if qt in ("cid", "identifier") and str(query).isdigit():
            return await self._search_by_cid(int(query))
        elif qt == "cid":
            return await self._search_by_cid(int(query))
        elif qt == "smiles":
            return await self._search_by_smiles(query)
        elif qt == "formula":
            return await self._search_by_formula(query)
        elif qt == "inchi":
            return await self._search_by_inchi(query)
        else:
            # Default: name search
            return await self._search_by_name(query)

    async def get_by_id(self, identifier: str) -> Compound | None:
        """Get a compound by CID.

        Args:
            identifier: PubChem CID.

        Returns:
            Compound object, or None if not found.
        """
        try:
            results = await self._search_by_cid(int(identifier))
            return results[0] if results else None
        except (NotFoundError, ValueError):
            return None

    async def get_properties(self, identifier: str) -> dict:
        """Get detailed properties for a compound by CID.

        Args:
            identifier: PubChem CID.

        Returns:
            Dictionary of property names to values.
        """
        path = f"compound/cid/{identifier}/property/{PROPERTY_LIST}/JSON"
        data = await self._http.get(path)
        props_list = data.get("PropertyTable", {}).get("Properties", [])
        if props_list:
            return props_list[0]
        return {}

    async def get_similarity(
        self,
        smiles: str,
        threshold: int = 90,
        max_results: int = 100,
    ) -> list[Compound]:
        """Find compounds similar to a given SMILES.

        Uses the PubChem similarity search API with Tanimoto threshold.

        Args:
            smiles: SMILES string to find similar compounds for.
            threshold: Tanimoto similarity threshold (0-100).
            max_results: Maximum number of results to return.

        Returns:
            List of similar Compound objects.
        """
        encoded = quote(smiles, safe="")
        path = f"compound/similarity/smiles/{encoded}/JSON"
        params = {"Threshold": threshold, "MaxRecords": max_results}

        try:
            data = await self._http.get(path, params=params, use_cache=False)
        except SourceError:
            return []

        # Handle ListKey async polling
        cids = await self._resolve_list_key_or_ids(data)
        if not cids:
            return []

        return await self._fetch_compounds_by_cids(cids[:max_results])

    async def get_substructure(
        self,
        smiles: str,
        max_results: int = 100,
    ) -> list[Compound]:
        """Find compounds containing a substructure.

        Args:
            smiles: SMILES substructure query.
            max_results: Maximum number of results.

        Returns:
            List of Compound objects containing the substructure.
        """
        encoded = quote(smiles, safe="")
        path = f"compound/substructure/smiles/{encoded}/JSON"
        params = {"MaxRecords": max_results}

        try:
            data = await self._http.get(path, params=params, use_cache=False)
        except SourceError:
            return []

        cids = await self._resolve_list_key_or_ids(data)
        if not cids:
            return []

        return await self._fetch_compounds_by_cids(cids[:max_results])

    async def get_bioassays(self, cid: int) -> list[dict[str, Any]]:
        """Get bioassay summary for a compound.

        Args:
            cid: PubChem Compound ID.

        Returns:
            List of bioassay summary dictionaries.
        """
        try:
            path = f"compound/cid/{cid}/assaysummary/JSON"
            data = await self._http.get(path)
            # Returns Table with AssaySummaries
            table = data.get("Table", {})
            columns = table.get("Columns", {}).get("Column", [])
            rows = table.get("Row", [])

            assays: list[dict[str, Any]] = []
            for row in rows:
                cells = row.get("Cell", [])
                if len(cells) == len(columns):
                    assays.append(dict(zip(columns, cells, strict=True)))
            return assays
        except (NotFoundError, SourceError):
            return []

    async def _search_by_name(self, name: str) -> list[Compound]:
        """Search by compound name."""
        encoded = quote(name, safe="")
        path = f"compound/name/{encoded}/property/{PROPERTY_LIST}/JSON"
        try:
            data = await self._http.get(path)
            compounds = self._parse_property_table(data)
            # Enrich with synonyms for the first result
            if compounds and compounds[0].cid:
                await self._enrich_synonyms(compounds)
            return compounds
        except NotFoundError as exc:
            raise NotFoundError(f"Compound not found: {name}", identifier=name) from exc

    async def _search_by_cid(self, cid: int) -> list[Compound]:
        """Search by CID."""
        path = f"compound/cid/{cid}/property/{PROPERTY_LIST}/JSON"
        try:
            data = await self._http.get(path)
            compounds = self._parse_property_table(data)
            if compounds:
                await self._enrich_synonyms(compounds)
            return compounds
        except NotFoundError as exc:
            raise NotFoundError(f"CID not found: {cid}", identifier=str(cid)) from exc

    async def _search_by_smiles(self, smiles: str) -> list[Compound]:
        """Search by SMILES (exact match)."""
        encoded = quote(smiles, safe="")
        path = f"compound/smiles/{encoded}/property/{PROPERTY_LIST}/JSON"
        try:
            data = await self._http.get(path)
            return self._parse_property_table(data)
        except NotFoundError as exc:
            raise NotFoundError(f"SMILES not found: {smiles}", identifier=smiles) from exc

    async def _search_by_formula(self, formula: str) -> list[Compound]:
        """Search by molecular formula (returns multiple results)."""
        path = f"compound/fastformula/{formula}/cids/JSON"
        try:
            data = await self._http.get(path, use_cache=False)
        except NotFoundError as exc:
            raise NotFoundError(f"Formula not found: {formula}", identifier=formula) from exc

        cids = await self._resolve_list_key_or_ids(data)
        if not cids:
            return []

        return await self._fetch_compounds_by_cids(cids[:50])

    async def _search_by_inchi(self, inchi: str) -> list[Compound]:
        """Search by InChI (using POST)."""
        path = f"compound/inchi/property/{PROPERTY_LIST}/JSON"
        try:
            data = await self._http.post_form(path, form_data={"inchi": inchi})
            return self._parse_property_table(data)
        except NotFoundError as exc:
            raise NotFoundError(f"InChI not found: {inchi}", identifier=inchi) from exc

    async def _fetch_compounds_by_cids(self, cids: list[int]) -> list[Compound]:
        """Fetch compound data for a list of CIDs."""
        if not cids:
            return []

        cid_str = ",".join(str(c) for c in cids)
        path = f"compound/cid/{cid_str}/property/{PROPERTY_LIST}/JSON"
        try:
            data = await self._http.get(path)
            return self._parse_property_table(data)
        except (NotFoundError, SourceError):
            return []

    async def _enrich_synonyms(self, compounds: list[Compound]) -> None:
        """Enrich compounds with synonyms from PubChem."""
        for compound in compounds[:5]:  # Limit to first 5 to avoid too many requests
            if compound.cid:
                try:
                    path = f"compound/cid/{compound.cid}/synonyms/JSON"
                    data = await self._http.get(path)
                    info_list = data.get("InformationList", {}).get("Information", [])
                    if info_list:
                        synonyms = info_list[0].get("Synonym", [])
                        compound.synonyms = synonyms[:20]  # Limit synonyms
                        if not compound.name and synonyms:
                            compound.name = synonyms[0]
                except (NotFoundError, SourceError):
                    pass

    async def _resolve_list_key_or_ids(self, data: Any) -> list[int]:
        """Resolve a PubChem response to a list of CIDs.

        Handles both direct CID lists and ListKey async polling.

        Args:
            data: Response data from PubChem.

        Returns:
            List of CIDs.
        """
        if not isinstance(data, dict):
            return []

        # Direct CID list
        if "IdentifierList" in data:
            return data["IdentifierList"].get("CID", [])

        # ListKey - async operation
        waiting = data.get("Waiting", {})
        if waiting:
            list_key = waiting.get("ListKey")
            if list_key:
                return await self._poll_list_key(list_key)

        return []

    async def _poll_list_key(
        self,
        list_key: str,
        max_wait: float = 120.0,
        poll_interval: float = 2.0,
    ) -> list[int]:
        """Poll a PubChem ListKey until results are ready.

        Args:
            list_key: The ListKey from the initial search request.
            max_wait: Maximum time to wait in seconds.
            poll_interval: Time between poll requests in seconds.

        Returns:
            List of CIDs from the completed search.

        Raises:
            TimeoutError: If the search does not complete within max_wait.
        """
        start_time = time.monotonic()
        path = f"compound/listkey/{list_key}/cids/JSON"

        while (time.monotonic() - start_time) < max_wait:
            try:
                result = await self._http.get(path, use_cache=False)
                if isinstance(result, dict):
                    waiting = result.get("Waiting", {})
                    if waiting:
                        await asyncio.sleep(poll_interval)
                        continue

                    id_list = result.get("IdentifierList", {})
                    cids: list[int] = id_list.get("CID", [])
                    if cids:
                        return cids

                    if "CID" in result:
                        return result["CID"]

                await asyncio.sleep(poll_interval)
            except (SourceError, NotFoundError):
                await asyncio.sleep(poll_interval)

        raise TimeoutError(
            f"PubChem async search timed out after {max_wait}s (ListKey: {list_key})"
        )

    def _parse_property_table(self, data: Any) -> list[Compound]:
        """Parse PubChem PropertyTable response into Compound objects.

        Args:
            data: Parsed JSON response from PubChem.

        Returns:
            List of Compound objects.
        """
        if not isinstance(data, dict):
            return []

        table = data.get("PropertyTable", {})
        props_list = table.get("Properties", [])

        compounds: list[Compound] = []
        for prop in props_list:
            compound = self._prop_dict_to_compound(prop)
            compounds.append(compound)

        return compounds

    def _prop_dict_to_compound(self, prop: dict[str, Any]) -> Compound:
        """Convert a PubChem property dict to a Compound.

        Args:
            prop: Property dictionary from PubChem API.

        Returns:
            Compound object.
        """
        properties = CompoundProperties(
            molecular_weight=prop.get("MolecularWeight"),
            exact_mass=prop.get("ExactMass"),
            xlogp=prop.get("XLogP"),
            tpsa=prop.get("TPSA"),
            hbd_count=prop.get("HBondDonorCount"),
            hba_count=prop.get("HBondAcceptorCount"),
            rotatable_bonds=prop.get("RotatableBondCount"),
            heavy_atom_count=prop.get("HeavyAtomCount"),
            complexity=prop.get("Complexity"),
        )

        smiles = (
            prop.get("CanonicalSMILES")
            or prop.get("SMILES")
            or prop.get("ConnectivitySMILES")
            or prop.get("IsomericSMILES")
            or ""
        )
        name = prop.get("IUPACName")

        return Compound(
            cid=prop.get("CID"),
            smiles=smiles,
            inchi=prop.get("InChI"),
            inchikey=prop.get("InChIKey"),
            name=name,
            formula=prop.get("MolecularFormula"),
            properties=properties,
            sources=["pubchem"],
        )

    def is_available(self) -> bool:
        """Check if PubChem is available (sync version checks URL format only)."""
        return True  # Async availability check not practical in sync context
