"""Tests for SureChEMBL patent chemistry source adapter."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.models.patent import Patent
from chemfuse.sources.surechembl import SureChEMBLAdapter

SURECHEMBL_BASE = "https://www.ebi.ac.uk/surechembl"
SURECHEMBL_API = f"{SURECHEMBL_BASE}/api/search"

# ---------------------------------------------------------------------------
# Mock response fixtures
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

SURECHEMBL_PATENT_RESPONSE = {
    "response": {
        "docs": [
            {
                "patent_id": "US-4568679-A",
                "title": "Coated aspirin tablet",
                "filing_date": "1984-03-15",
                "assignee": "Bayer AG",
                "jurisdiction": "US",
            }
        ]
    }
}

SURECHEMBL_MULTIPLE_PATENTS_RESPONSE = {
    "response": {
        "docs": [
            {
                "patent_id": "US-4568679-A",
                "title": "Coated aspirin tablet",
                "filing_date": "1984-03-15",
                "assignee": "Bayer AG",
                "jurisdiction": "US",
            },
            {
                "patent_id": "EP-0123456-B1",
                "title": "Aspirin formulation",
                "filing_date": "1990-06-20",
                "assignee": "Pharma Corp",
                "jurisdiction": "EP",
            },
        ]
    }
}

SURECHEMBL_EMPTY_RESPONSE = {
    "response": {
        "docs": []
    }
}

SURECHEMBL_FLAT_LIST_RESPONSE = [
    {
        "patent_id": "US-9999999-A",
        "title": "Novel aspirin derivative",
        "filing_date": "2020-01-01",
        "assignee": "University X",
        "jurisdiction": "US",
    }
]


class TestSureChEMBLSearchBySmiles:
    """Tests for SureChEMBLAdapter.search_by_smiles."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_returns_patent_objects(self) -> None:
        """search_by_smiles returns Patent objects."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_PATENT_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert len(patents) == 1
        assert all(isinstance(p, Patent) for p in patents)

    @pytest.mark.asyncio
    @respx.mock
    async def test_patent_fields_populated(self) -> None:
        """search_by_smiles correctly populates Patent fields."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_PATENT_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles(ASPIRIN_SMILES)
        p = patents[0]
        assert p.patent_id == "US-4568679-A"
        assert p.title == "Coated aspirin tablet"
        assert p.filing_date == "1984-03-15"
        assert p.assignee == "Bayer AG"
        assert p.jurisdiction == "US"

    @pytest.mark.asyncio
    @respx.mock
    async def test_returns_multiple_patents(self) -> None:
        """search_by_smiles returns multiple patents when available."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_MULTIPLE_PATENTS_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert len(patents) == 2

    @pytest.mark.asyncio
    @respx.mock
    async def test_empty_result_when_no_patents(self) -> None:
        """search_by_smiles returns empty list when no patents found."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_EMPTY_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles("CC")
        assert patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_api_error_returns_empty(self) -> None:
        """search_by_smiles returns empty list on API error."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(503, text="Service Unavailable")
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_source_is_surechembl(self) -> None:
        """All returned patents have source='surechembl'."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_PATENT_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert all(p.source == "surechembl" for p in patents)


class TestSureChEMBLSearchByInchikey:
    """Tests for SureChEMBLAdapter.search_by_inchikey."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_returns_patent_objects(self) -> None:
        """search_by_inchikey returns Patent objects."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_PATENT_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_inchikey(ASPIRIN_INCHIKEY)
        assert len(patents) == 1
        assert isinstance(patents[0], Patent)

    @pytest.mark.asyncio
    @respx.mock
    async def test_empty_result_when_not_found(self) -> None:
        """search_by_inchikey returns empty list when compound not in patents."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_EMPTY_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_inchikey("XXXXXXXXXXXXXXXX-YYYYYYYYYY-N")
        assert patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_api_error_returns_empty(self) -> None:
        """search_by_inchikey returns empty list on API error."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(500, text="Error")
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_inchikey(ASPIRIN_INCHIKEY)
        assert patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_patent_fields_match_smiles_response(self) -> None:
        """search_by_inchikey parses the same fields as SMILES search."""
        respx.get(SURECHEMBL_API).mock(
            return_value=Response(200, json=SURECHEMBL_PATENT_RESPONSE)
        )
        adapter = SureChEMBLAdapter()
        patents = await adapter.search_by_inchikey(ASPIRIN_INCHIKEY)
        p = patents[0]
        assert p.patent_id == "US-4568679-A"
        assert p.title == "Coated aspirin tablet"


class TestSureChEMBLResponseParsing:
    """Tests for SureChEMBL response parsing edge cases."""

    def test_parse_flat_list_response(self) -> None:
        """_parse_response handles flat list input."""
        adapter = SureChEMBLAdapter()
        patents = adapter._parse_response(SURECHEMBL_FLAT_LIST_RESPONSE)
        assert len(patents) == 1
        assert patents[0].patent_id == "US-9999999-A"

    def test_parse_response_none_returns_empty(self) -> None:
        """_parse_response returns empty list for None input."""
        adapter = SureChEMBLAdapter()
        result = adapter._parse_response(None)
        assert result == []

    def test_parse_response_empty_dict_returns_empty(self) -> None:
        """_parse_response returns empty list for empty dict."""
        adapter = SureChEMBLAdapter()
        result = adapter._parse_response({})
        assert result == []

    def test_parse_response_with_patents_key(self) -> None:
        """_parse_response handles 'patents' key in response dict."""
        adapter = SureChEMBLAdapter()
        data = {"patents": [{"patent_id": "US-111", "title": "Test"}]}
        patents = adapter._parse_response(data)
        assert len(patents) == 1
        assert patents[0].patent_id == "US-111"

    def test_parse_response_with_results_key(self) -> None:
        """_parse_response handles 'results' key in response dict."""
        adapter = SureChEMBLAdapter()
        data = {"results": [{"patent_id": "EP-222", "title": "Another test"}]}
        patents = adapter._parse_response(data)
        assert len(patents) == 1
        assert patents[0].patent_id == "EP-222"

    def test_parse_patent_doc_missing_patent_id_returns_none(self) -> None:
        """_parse_patent_doc returns None when patent_id is missing."""
        adapter = SureChEMBLAdapter()
        result = adapter._parse_patent_doc({"title": "No ID here"})
        assert result is None

    def test_parse_patent_doc_with_schembl_id(self) -> None:
        """_parse_patent_doc uses schembl_id field as fallback for patent_id."""
        adapter = SureChEMBLAdapter()
        doc = {
            "schembl_id": "SCHEMBL19",
            "patent_id": "US-4568679-A",
            "title": "Coated aspirin tablet",
        }
        patent = adapter._parse_patent_doc(doc)
        assert patent is not None
        assert patent.patent_id == "US-4568679-A"

    def test_parse_docs_skips_non_dict_items(self) -> None:
        """_parse_docs skips non-dict items in the list."""
        adapter = SureChEMBLAdapter()
        docs = [
            "not_a_dict",
            {"patent_id": "US-111"},
            42,
        ]
        patents = adapter._parse_docs(docs)
        assert len(patents) == 1
        assert patents[0].patent_id == "US-111"

    def test_parse_response_with_malformed_doc_skipped(self) -> None:
        """Docs with invalid data are skipped gracefully."""
        adapter = SureChEMBLAdapter()
        data = {
            "response": {
                "docs": [
                    {},  # missing patent_id
                    {"patent_id": "US-VALID"},
                ]
            }
        }
        patents = adapter._parse_response(data)
        assert len(patents) == 1
        assert patents[0].patent_id == "US-VALID"


class TestSureChEMBLAdapterMisc:
    """Miscellaneous SureChEMBL adapter tests."""

    def test_adapter_name(self) -> None:
        """Adapter has correct name."""
        adapter = SureChEMBLAdapter()
        assert adapter.name == "surechembl"

    def test_is_available(self) -> None:
        """is_available always returns True."""
        adapter = SureChEMBLAdapter()
        assert adapter.is_available() is True

    @pytest.mark.asyncio
    async def test_search_returns_empty(self) -> None:
        """search always returns empty list."""
        adapter = SureChEMBLAdapter()
        result = await adapter.search("aspirin")
        assert result == []

    @pytest.mark.asyncio
    async def test_get_by_id_returns_none(self) -> None:
        """get_by_id always returns None."""
        adapter = SureChEMBLAdapter()
        result = await adapter.get_by_id("SCHEMBL19")
        assert result is None

    @pytest.mark.asyncio
    async def test_get_properties_returns_empty(self) -> None:
        """get_properties always returns empty dict."""
        adapter = SureChEMBLAdapter()
        result = await adapter.get_properties("SCHEMBL19")
        assert result == {}

    def test_rate_limit_attribute(self) -> None:
        """Adapter has rate_limit attribute."""
        adapter = SureChEMBLAdapter()
        assert adapter.rate_limit > 0
