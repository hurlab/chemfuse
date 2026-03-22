"""UniChem cross-reference adapter."""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.models.compound import Compound
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

UNICHEM_BASE_URL = "https://www.ebi.ac.uk/unichem/rest"
UNICHEM_RATE_LIMIT = 5.0  # no hard limit; use a conservative rate

logger = logging.getLogger(__name__)

# UniChem source IDs for major databases
SOURCE_IDS: dict[str, int] = {
    "chembl": 1,
    "drugbank": 2,
    "kegg": 6,
    "chebi": 7,
    "pubchem": 22,
}

# Reverse mapping: ID -> database name
SOURCE_NAMES: dict[int, str] = {v: k for k, v in SOURCE_IDS.items()}


class UniChemAdapter(SourceAdapter):
    """UniChem REST API adapter for cross-database identifier mapping.

    Maps compound identifiers between major chemical databases using
    InChIKey-based or source-specific identifier lookup.

    No API key required.
    Base URL: https://www.ebi.ac.uk/unichem/rest/
    """

    name = "unichem"
    base_url = UNICHEM_BASE_URL
    rate_limit = UNICHEM_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the UniChem adapter.

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
        """Search is not supported for UniChem; returns empty list.

        UniChem is a cross-reference tool, not a compound search engine.
        Use map_identifiers() or cross_reference() instead.
        """
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        """Get by ID is not applicable for UniChem.

        Returns:
            None always; use cross_reference() instead.
        """
        return None

    async def get_properties(self, identifier: str) -> dict:
        """Get properties is not supported for UniChem.

        Returns:
            Empty dict.
        """
        return {}

    def is_available(self) -> bool:
        """Always returns True."""
        return True

    async def map_identifiers(
        self,
        identifier: str,
        source_type: str,
    ) -> dict[str, str]:
        """Map an identifier from one database to all other databases.

        Args:
            identifier: The compound identifier to map.
            source_type: Source database name (e.g., "pubchem", "chembl").

        Returns:
            Dictionary mapping database names to identifiers.
            Returns empty dict if no mappings found.
        """
        src_id = SOURCE_IDS.get(source_type.lower())
        if src_id is None:
            logger.warning("Unknown UniChem source type: %s", source_type)
            return {}

        try:
            data = await self._http.get(
                f"/src_compound_id/{identifier}/src_id/{src_id}",
            )
        except Exception as exc:
            logger.debug("UniChem map_identifiers failed: %s", exc)
            return {}

        return self._parse_mappings(data)

    async def cross_reference(self, inchikey: str) -> dict[str, str]:
        """Get all database identifiers for a compound by InChIKey.

        Args:
            inchikey: Standard InChIKey (27 characters).

        Returns:
            Dictionary mapping database names to identifiers.
            Returns empty dict if no mappings found.
        """
        try:
            data = await self._http.get(f"/inchikey/{inchikey}")
        except Exception as exc:
            logger.debug("UniChem cross_reference failed for %s: %s", inchikey, exc)
            return {}

        return self._parse_mappings(data)

    async def batch_map(
        self,
        identifiers: list[str],
        source_type: str,
    ) -> dict[str, dict[str, str]]:
        """Map multiple identifiers in a single batch operation.

        Args:
            identifiers: List of compound identifiers.
            source_type: Source database name.

        Returns:
            Dictionary mapping each input identifier to its cross-reference dict.
        """
        results: dict[str, dict[str, str]] = {}
        for ident in identifiers:
            results[ident] = await self.map_identifiers(ident, source_type)
        return results

    # --- Private helpers ---

    def _parse_mappings(self, data: Any) -> dict[str, str]:
        """Parse UniChem mapping response into a database-name -> ID dict.

        Args:
            data: Raw response from UniChem API (list of dicts).

        Returns:
            Mapping from database names to identifier strings.
        """
        if not isinstance(data, list):
            return {}

        result: dict[str, str] = {}
        for entry in data:
            if not isinstance(entry, dict):
                continue
            src_id_raw = entry.get("src_id")
            compound_id = entry.get("src_compound_id") or entry.get("compound_id")
            if src_id_raw is None or compound_id is None:
                continue
            try:
                src_id = int(src_id_raw)
            except (TypeError, ValueError):
                continue
            db_name = SOURCE_NAMES.get(src_id)
            if db_name:
                result[db_name] = str(compound_id)
        return result


# backward compatibility alias
UniChemSource = UniChemAdapter

__all__ = ["UniChemAdapter", "UniChemSource", "SOURCE_IDS", "SOURCE_NAMES"]
