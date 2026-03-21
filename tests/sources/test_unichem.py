"""Tests for UniChem cross-reference adapter."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.sources.unichem import SOURCE_IDS, SOURCE_NAMES, UniChemAdapter

UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/rest"

# Aspirin InChIKey
ASPIRIN_INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

# Mock UniChem response for aspirin InChIKey lookup
ASPIRIN_UNICHEM_RESPONSE = [
    {"src_id": "1", "src_compound_id": "CHEMBL25"},      # ChEMBL
    {"src_id": "22", "src_compound_id": "2244"},           # PubChem
    {"src_id": "2", "src_compound_id": "DB00945"},         # DrugBank
    {"src_id": "6", "src_compound_id": "D00109"},          # KEGG
    {"src_id": "7", "src_compound_id": "CHEBI:15365"},     # ChEBI
]

# Response from mapping PubChem CID 2244
PUBCHEM_MAP_RESPONSE = [
    {"src_id": "1", "src_compound_id": "CHEMBL25"},
    {"src_id": "2", "src_compound_id": "DB00945"},
]


class TestUniChemSourceConstants:
    """Tests for UniChem source ID constants."""

    def test_source_ids_contain_required_databases(self) -> None:
        """SOURCE_IDS contains required database mappings."""
        assert SOURCE_IDS["pubchem"] == 22
        assert SOURCE_IDS["chembl"] == 1
        assert SOURCE_IDS["drugbank"] == 2
        assert SOURCE_IDS["kegg"] == 6
        assert SOURCE_IDS["chebi"] == 7

    def test_source_names_is_inverse_of_source_ids(self) -> None:
        """SOURCE_NAMES is the inverse mapping of SOURCE_IDS."""
        for name, src_id in SOURCE_IDS.items():
            assert SOURCE_NAMES[src_id] == name


class TestUniChemCrossReference:
    """Tests for UniChem InChIKey-based cross-reference lookup."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_cross_reference_by_inchikey(self) -> None:
        """cross_reference maps InChIKey to all database IDs."""
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=ASPIRIN_UNICHEM_RESPONSE)
        )
        adapter = UniChemAdapter()
        result = await adapter.cross_reference(ASPIRIN_INCHIKEY)
        assert result["chembl"] == "CHEMBL25"
        assert result["pubchem"] == "2244"
        assert result["drugbank"] == "DB00945"
        assert result["kegg"] == "D00109"
        assert result["chebi"] == "CHEBI:15365"

    @pytest.mark.asyncio
    @respx.mock
    async def test_cross_reference_empty_on_not_found(self) -> None:
        """cross_reference returns empty dict when no mappings found."""
        respx.get(f"{UNICHEM_BASE}/inchikey/UNKNOWNKEY-UHFFFAOYSA-N").mock(
            return_value=Response(404, json={"error": "not found"})
        )
        adapter = UniChemAdapter()
        result = await adapter.cross_reference("UNKNOWNKEY-UHFFFAOYSA-N")
        assert result == {}

    @pytest.mark.asyncio
    @respx.mock
    async def test_cross_reference_empty_list_response(self) -> None:
        """cross_reference returns empty dict on empty list response."""
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=[])
        )
        adapter = UniChemAdapter()
        result = await adapter.cross_reference(ASPIRIN_INCHIKEY)
        assert result == {}

    @pytest.mark.asyncio
    @respx.mock
    async def test_cross_reference_ignores_unknown_source_ids(self) -> None:
        """cross_reference skips entries with unknown source IDs."""
        data = [
            {"src_id": "99", "src_compound_id": "UNKNOWN123"},  # Unknown
            {"src_id": "1", "src_compound_id": "CHEMBL25"},     # Known
        ]
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=data)
        )
        adapter = UniChemAdapter()
        result = await adapter.cross_reference(ASPIRIN_INCHIKEY)
        assert "chembl" in result
        assert len(result) == 1  # Only chembl, unknown ID ignored


class TestUniChemMapIdentifiers:
    """Tests for UniChem map_identifiers method."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_pubchem_cid_to_chembl(self) -> None:
        """Map PubChem CID 2244 to ChEMBL ID."""
        respx.get(f"{UNICHEM_BASE}/src_compound_id/2244/src_id/22").mock(
            return_value=Response(200, json=PUBCHEM_MAP_RESPONSE)
        )
        adapter = UniChemAdapter()
        result = await adapter.map_identifiers("2244", "pubchem")
        assert result["chembl"] == "CHEMBL25"
        assert result["drugbank"] == "DB00945"

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_chembl_id_to_pubchem(self) -> None:
        """Map ChEMBL ID to PubChem CID."""
        respx.get(f"{UNICHEM_BASE}/src_compound_id/CHEMBL25/src_id/1").mock(
            return_value=Response(200, json=[
                {"src_id": "22", "src_compound_id": "2244"},
            ])
        )
        adapter = UniChemAdapter()
        result = await adapter.map_identifiers("CHEMBL25", "chembl")
        assert result["pubchem"] == "2244"

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_unknown_source_returns_empty(self) -> None:
        """map_identifiers returns empty dict for unknown source type."""
        adapter = UniChemAdapter()
        result = await adapter.map_identifiers("12345", "unknown_source")
        assert result == {}

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_not_found_returns_empty(self) -> None:
        """map_identifiers returns empty dict on API error."""
        respx.get(f"{UNICHEM_BASE}/src_compound_id/999/src_id/22").mock(
            return_value=Response(404, json={"error": "not found"})
        )
        adapter = UniChemAdapter()
        result = await adapter.map_identifiers("999", "pubchem")
        assert result == {}


class TestUniChemBatchMap:
    """Tests for UniChem batch_map method."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_batch_map_multiple_identifiers(self) -> None:
        """batch_map maps each identifier independently."""
        respx.get(f"{UNICHEM_BASE}/src_compound_id/2244/src_id/22").mock(
            return_value=Response(200, json=[{"src_id": "1", "src_compound_id": "CHEMBL25"}])
        )
        respx.get(f"{UNICHEM_BASE}/src_compound_id/3672/src_id/22").mock(
            return_value=Response(200, json=[{"src_id": "1", "src_compound_id": "CHEMBL521"}])
        )
        adapter = UniChemAdapter()
        results = await adapter.batch_map(["2244", "3672"], "pubchem")
        assert results["2244"]["chembl"] == "CHEMBL25"
        assert results["3672"]["chembl"] == "CHEMBL521"

    @pytest.mark.asyncio
    @respx.mock
    async def test_batch_map_empty_list(self) -> None:
        """batch_map with empty list returns empty dict."""
        adapter = UniChemAdapter()
        results = await adapter.batch_map([], "pubchem")
        assert results == {}


class TestUniChemAdapterMisc:
    """Miscellaneous UniChem adapter tests."""

    def test_adapter_name(self) -> None:
        """Adapter has correct name."""
        adapter = UniChemAdapter()
        assert adapter.name == "unichem"

    def test_is_available(self) -> None:
        """is_available always returns True."""
        adapter = UniChemAdapter()
        assert adapter.is_available() is True

    @pytest.mark.asyncio
    async def test_search_returns_empty(self) -> None:
        """search always returns empty list."""
        adapter = UniChemAdapter()
        result = await adapter.search("aspirin")
        assert result == []

    @pytest.mark.asyncio
    async def test_get_by_id_returns_none(self) -> None:
        """get_by_id always returns None."""
        adapter = UniChemAdapter()
        result = await adapter.get_by_id("CHEMBL25")
        assert result is None

    @pytest.mark.asyncio
    async def test_get_properties_returns_empty(self) -> None:
        """get_properties always returns empty dict."""
        adapter = UniChemAdapter()
        result = await adapter.get_properties("CHEMBL25")
        assert result == {}

    @pytest.mark.asyncio
    @respx.mock
    async def test_malformed_response_returns_empty(self) -> None:
        """Malformed (non-list) response returns empty dict."""
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json={"unexpected": "format"})
        )
        adapter = UniChemAdapter()
        result = await adapter.cross_reference(ASPIRIN_INCHIKEY)
        assert result == {}
