"""SureChEMBL patent chemistry source adapter."""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.models.compound import Compound
from chemfuse.models.patent import Patent
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

SURECHEMBL_BASE_URL = "https://www.ebi.ac.uk/surechembl"
SURECHEMBL_RATE_LIMIT = 2.0  # conservative; no hard limit documented

logger = logging.getLogger(__name__)


class SureChEMBLAdapter(SourceAdapter):
    """SureChEMBL patent chemistry source adapter.

    Note: This adapter is experimental. SureChEMBL does not provide a public
    REST API, so functionality is limited.

    Searches patent literature for compounds by SMILES or InChIKey.
    No API key required (CC BY 4.0 license, public access).

    API base: https://www.ebi.ac.uk/surechembl
    """

    name = "surechembl"
    base_url = SURECHEMBL_BASE_URL
    rate_limit = SURECHEMBL_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the SureChEMBL adapter.

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
        query_type: str = "smiles",
    ) -> list[Compound]:
        """Search is a pass-through for SureChEMBL.

        SureChEMBL is used for patent data enrichment. Returns empty list.
        Use search_by_smiles() or search_by_inchikey() directly.
        """
        logger.warning("SureChEMBL adapter is experimental — limited functionality")
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        """Not implemented for SureChEMBL.

        Returns:
            None always.
        """
        logger.warning("SureChEMBL adapter is experimental — limited functionality")
        return None

    async def get_properties(self, identifier: str) -> dict:
        """Not implemented for SureChEMBL.

        Returns:
            Empty dict.
        """
        logger.warning("SureChEMBL adapter is experimental — limited functionality")
        return {}

    def is_available(self) -> bool:
        """Always returns True."""
        return True

    async def search_by_smiles(
        self,
        smiles: str,
    ) -> list[Patent]:
        """Search for patents containing a compound by SMILES.

        Args:
            smiles: SMILES string for the query compound.

        Returns:
            List of Patent objects, empty list if none found.
        """
        try:
            data = await self._http.get(
                "/api/search",
                params={"smiles": smiles, "limit": 100},
            )
        except Exception as exc:
            logger.debug("SureChEMBL search_by_smiles failed: %s", exc)
            return []

        return self._parse_response(data)

    async def search_by_inchikey(
        self,
        inchikey: str,
    ) -> list[Patent]:
        """Search for patents containing a compound by InChIKey.

        Args:
            inchikey: Standard InChIKey for the compound.

        Returns:
            List of Patent objects, empty list if none found.
        """
        try:
            data = await self._http.get(
                "/api/search",
                params={"inchikey": inchikey, "limit": 100},
            )
        except Exception as exc:
            logger.debug("SureChEMBL search_by_inchikey failed: %s", exc)
            return []

        return self._parse_response(data)

    # --- Private helpers ---

    def _parse_response(self, data: Any) -> list[Patent]:
        """Parse SureChEMBL API response into Patent objects.

        Handles the standard SureChEMBL JSON response format:
        {"response": {"docs": [...patent records...]}}

        Also handles flat lists and direct patent arrays.

        Args:
            data: Raw API response.

        Returns:
            List of Patent objects.
        """
        if data is None:
            return []

        # Standard SureChEMBL format: {"response": {"docs": [...]}}
        if isinstance(data, dict):
            response = data.get("response")
            if isinstance(response, dict):
                docs = response.get("docs") or []
                return self._parse_docs(docs)
            # Direct list under a 'patents' or 'results' key
            for key in ("patents", "results", "docs"):
                items = data.get(key)
                if isinstance(items, list):
                    return self._parse_docs(items)
            return []

        if isinstance(data, list):
            return self._parse_docs(data)

        return []

    def _parse_docs(self, docs: list) -> list[Patent]:
        """Parse a list of patent document dicts into Patent objects.

        Args:
            docs: List of raw patent document dicts.

        Returns:
            List of Patent objects.
        """
        patents: list[Patent] = []
        for doc in docs:
            if not isinstance(doc, dict):
                continue
            patent = self._parse_patent_doc(doc)
            if patent is not None:
                patents.append(patent)
        return patents

    def _parse_patent_doc(self, doc: dict) -> Patent | None:
        """Convert a single patent document dict to a Patent object.

        Args:
            doc: Raw patent document dict from SureChEMBL.

        Returns:
            Patent object or None if the patent_id is missing.
        """
        patent_id = (
            doc.get("patent_id")
            or doc.get("patentId")
            or doc.get("id")
            or doc.get("schembl_id")
            or ""
        )
        if not patent_id or not str(patent_id).strip():
            return None

        title = (
            doc.get("title")
            or doc.get("patent_title")
            or None
        )
        filing_date = (
            doc.get("filing_date")
            or doc.get("filingDate")
            or doc.get("date")
            or None
        )
        assignee = (
            doc.get("assignee")
            or doc.get("patent_assignee")
            or None
        )
        jurisdiction = (
            doc.get("jurisdiction")
            or doc.get("country")
            or None
        )

        try:
            return Patent(
                patent_id=str(patent_id).strip(),
                title=str(title) if title else None,
                filing_date=str(filing_date) if filing_date else None,
                assignee=str(assignee) if assignee else None,
                jurisdiction=str(jurisdiction) if jurisdiction else None,
                source=self.name,
            )
        except Exception as exc:
            logger.debug("Failed to parse patent doc: %s - %s", doc, exc)
            return None


__all__ = ["SureChEMBLAdapter"]
