"""BindingDB REST API source adapter."""

from __future__ import annotations

import logging
import re
from typing import Any

import defusedxml.ElementTree as ET

from chemfuse.models.bioactivity import BindingMeasurement
from chemfuse.models.compound import Compound
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

BINDINGDB_BASE_URL = "https://www.bindingdb.org/axis2/services/BDBService"
BINDINGDB_RATE_LIMIT = 1.0  # requests per second

logger = logging.getLogger(__name__)

# Pattern to extract numeric value with optional inequality prefix
_INEQUALITY_RE = re.compile(r"^([<>]=?|~)?\s*([\d.]+)\s*$")


def _parse_affinity(raw: Any) -> tuple[float | None, str | None]:
    """Parse an affinity value that may include an inequality prefix.

    Args:
        raw: Raw value string or number (e.g., ">10000", "1.67", "~100").

    Returns:
        Tuple of (numeric_value, relation_string).
        relation_string is None for plain numeric values.
    """
    if raw is None:
        return None, None
    text = str(raw).strip()
    if not text or text.lower() in ("", "na", "n/a", "none", "-"):
        return None, None
    match = _INEQUALITY_RE.match(text)
    if match:
        relation = match.group(1)
        value = float(match.group(2))
        return value, relation or None
    # Fallback: try plain float
    try:
        return float(text), None
    except ValueError:
        return None, None


class BindingDBAdapter(SourceAdapter):
    """BindingDB REST API adapter for binding affinity measurements.

    Provides search by SMILES string and by UniProt target ID.
    Returns Ki, Kd, IC50, and EC50 measurements.

    Rate limit: 1 req/sec (no API key required).
    API base: https://www.bindingdb.org/axis2/services/BDBService/
    """

    name = "bindingdb"
    base_url = BINDINGDB_BASE_URL
    rate_limit = BINDINGDB_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the BindingDB adapter.

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
        """Search is minimal for BindingDB; returns empty list.

        BindingDB is used for binding measurements. Use search_by_smiles()
        or search_by_target() for binding data retrieval.
        """
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        """Not implemented for BindingDB.

        Returns:
            None always.
        """
        return None

    async def get_properties(self, identifier: str) -> dict:
        """Not implemented for BindingDB.

        Returns:
            Empty dict.
        """
        return {}

    def is_available(self) -> bool:
        """Always returns True."""
        return True

    async def search_by_smiles(
        self,
        smiles: str,
    ) -> list[BindingMeasurement]:
        """Find binding measurements for a compound by SMILES.

        Args:
            smiles: SMILES string for the query compound.

        Returns:
            List of BindingMeasurement objects, empty list if none found.
        """
        try:
            data = await self._http.get(
                "/getLigandsBySmiles",
                params={"smiles": smiles, "response": "application/json"},
            )
        except Exception as exc:
            logger.debug("BindingDB search_by_smiles failed: %s", exc)
            return []

        return self._parse_response(data)

    async def search_by_target(
        self,
        uniprot_id: str,
    ) -> list[BindingMeasurement]:
        """Find binding measurements for compounds targeting a specific protein.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P23219").

        Returns:
            List of BindingMeasurement objects, empty list if none found.
        """
        try:
            data = await self._http.get(
                "/getLigandsByUniProt",
                params={"uniprot": uniprot_id, "response": "application/json"},
            )
        except Exception as exc:
            logger.debug("BindingDB search_by_target failed for %s: %s", uniprot_id, exc)
            return []

        return self._parse_response(data)

    # --- Private helpers ---

    def _parse_response(self, data: Any) -> list[BindingMeasurement]:
        """Parse BindingDB response (JSON or XML string) into measurements.

        Args:
            data: Raw response data (dict, list, or XML string).

        Returns:
            List of BindingMeasurement objects.
        """
        if isinstance(data, str):
            # Try XML parsing
            return self._parse_xml(data)
        if isinstance(data, dict):
            return self._parse_json(data)
        if isinstance(data, list):
            return self._parse_json_list(data)
        return []

    def _parse_json(self, data: dict) -> list[BindingMeasurement]:
        """Parse JSON dict response from BindingDB."""
        # BindingDB may return data under various keys
        measurements_raw = (
            data.get("affinities")
            or data.get("ligands")
            or data.get("results")
            or data.get("getLigandsBySmilesResponse")
            or data.get("getLigandsByUniProtResponse")
            or []
        )
        if isinstance(measurements_raw, dict):
            # Unwrap nested structure
            measurements_raw = (
                measurements_raw.get("affinities")
                or measurements_raw.get("ligands")
                or []
            )
        if not isinstance(measurements_raw, list):
            return []
        return self._parse_json_list(measurements_raw)

    def _parse_json_list(self, items: list) -> list[BindingMeasurement]:
        """Parse a list of measurement dicts from BindingDB."""
        result = []
        for item in items:
            if not isinstance(item, dict):
                continue
            measurement = self._parse_measurement_dict(item)
            if measurement is not None:
                result.append(measurement)
        return result

    def _parse_measurement_dict(self, item: dict) -> BindingMeasurement | None:
        """Convert a single measurement dict to a BindingMeasurement.

        Args:
            item: Raw measurement dict.

        Returns:
            BindingMeasurement or None if data is insufficient.
        """
        target_name = (
            item.get("target")
            or item.get("target_name")
            or item.get("name")
            or ""
        )
        uniprot = item.get("uniprot") or item.get("target_uniprot") or item.get("uniprot_id")
        reference = item.get("doi") or item.get("pmid") or item.get("source")

        ki_raw = item.get("ki") or item.get("Ki")
        kd_raw = item.get("kd") or item.get("Kd")
        ic50_raw = item.get("ic50") or item.get("IC50")
        ec50_raw = item.get("ec50") or item.get("EC50")

        ki, ki_rel = _parse_affinity(ki_raw)
        kd, kd_rel = _parse_affinity(kd_raw)
        ic50, ic50_rel = _parse_affinity(ic50_raw)
        ec50, ec50_rel = _parse_affinity(ec50_raw)

        # Skip records with no target name and no measurements
        if not target_name and all(v is None for v in [ki, kd, ic50, ec50]):
            return None

        return BindingMeasurement(
            target_name=target_name or "Unknown",
            target_uniprot=uniprot,
            ki=ki,
            kd=kd,
            ic50=ic50,
            ec50=ec50,
            ki_relation=ki_rel,
            kd_relation=kd_rel,
            ic50_relation=ic50_rel,
            ec50_relation=ec50_rel,
            source="bindingdb",
            reference=str(reference) if reference else None,
        )

    def _parse_xml(self, xml_text: str) -> list[BindingMeasurement]:
        """Parse XML response from BindingDB.

        Args:
            xml_text: XML string from BindingDB API.

        Returns:
            List of BindingMeasurement objects.
        """
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as exc:
            logger.debug("BindingDB XML parse error: %s", exc)
            return []

        result = []
        # Handle various XML structures
        # Look for affinities or ligands elements anywhere in the tree
        for affinity in root.iter("affinity"):
            target_el = affinity.find("target")
            target_name = target_el.text.strip() if target_el is not None and target_el.text else ""

            uniprot_el = affinity.find("uniprot")
            uniprot = uniprot_el.text.strip() if uniprot_el is not None and uniprot_el.text else None

            ki_el = affinity.find("ki")
            kd_el = affinity.find("kd")
            ic50_el = affinity.find("ic50")
            ec50_el = affinity.find("ec50")

            ki, ki_rel = _parse_affinity(ki_el.text if ki_el is not None else None)
            kd, kd_rel = _parse_affinity(kd_el.text if kd_el is not None else None)
            ic50, ic50_rel = _parse_affinity(ic50_el.text if ic50_el is not None else None)
            ec50, ec50_rel = _parse_affinity(ec50_el.text if ec50_el is not None else None)

            if not target_name and all(v is None for v in [ki, kd, ic50, ec50]):
                continue

            result.append(BindingMeasurement(
                target_name=target_name or "Unknown",
                target_uniprot=uniprot,
                ki=ki,
                kd=kd,
                ic50=ic50,
                ec50=ec50,
                ki_relation=ki_rel,
                kd_relation=kd_rel,
                ic50_relation=ic50_rel,
                ec50_relation=ec50_rel,
                source="bindingdb",
            ))

        return result


__all__ = ["BindingDBAdapter"]
