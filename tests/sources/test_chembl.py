"""Tests for ChEMBL source adapter."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.models.bioactivity import Bioactivity
from chemfuse.sources.chembl import ChEMBLAdapter

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

# --- Mock response fixtures ---

ASPIRIN_MOLECULE = {
    "molecule_chembl_id": "CHEMBL25",
    "pref_name": "ASPIRIN",
    "molecule_structures": {
        "canonical_smiles": "CC(=O)Oc1ccccc1C(O)=O",
        "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    },
    "molecule_properties": {
        "full_mwt": "180.16",
        "alogp": "1.31",
        "hbd": "1",
        "hba": "4",
        "psa": "63.60",
        "rtb": "3",
        "heavy_atoms": "13",
        "molecular_formula": "C9H8O4",
    },
}

ASPIRIN_SEARCH_RESPONSE = {
    "molecules": [ASPIRIN_MOLECULE],
    "page_meta": {
        "total_count": 1,
        "limit": 20,
        "offset": 0,
        "next": None,
        "previous": None,
    },
}

ASPIRIN_ACTIVITIES_PAGE1 = {
    "activities": [
        {
            "target_pref_name": "Cyclooxygenase-1",
            "target_chembl_id": "CHEMBL247",
            "standard_type": "IC50",
            "standard_value": "1670",
            "standard_units": "nM",
            "standard_relation": "=",
            "assay_type": "B",
            "document_chembl_id": "CHEMBL1234",
        },
        {
            "target_pref_name": "Cyclooxygenase-2",
            "target_chembl_id": "CHEMBL230",
            "standard_type": "IC50",
            "standard_value": "278000",
            "standard_units": "nM",
            "standard_relation": ">",
            "assay_type": "B",
            "document_chembl_id": "CHEMBL5678",
        },
    ],
    "page_meta": {
        "total_count": 3,
        "limit": 2,
        "offset": 0,
        "next": f"{CHEMBL_BASE}/activity?molecule_chembl_id=CHEMBL25&format=json&limit=2&offset=2",
        "previous": None,
    },
}

ASPIRIN_ACTIVITIES_PAGE2 = {
    "activities": [
        {
            "target_pref_name": "Prostaglandin E2 receptor",
            "target_chembl_id": "CHEMBL290",
            "standard_type": "Ki",
            "standard_value": "500",
            "standard_units": "nM",
            "standard_relation": "=",
            "assay_type": "B",
        },
    ],
    "page_meta": {
        "total_count": 3,
        "limit": 2,
        "offset": 2,
        "next": None,
        "previous": f"{CHEMBL_BASE}/activity?molecule_chembl_id=CHEMBL25&format=json&limit=2&offset=0",
    },
}

ASPIRIN_MECHANISMS = {
    "mechanisms": [
        {
            "target_chembl_id": "CHEMBL247",
            "action_type": "INHIBITOR",
            "mechanism_of_action": "Cyclooxygenase inhibitor",
            "mechanism_refs": [{"ref_type": "PubMed", "ref_id": "12345"}],
        }
    ]
}


class TestChEMBLSearch:
    """Tests for ChEMBL search functionality."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_name(self) -> None:
        """Search by name returns Compound objects."""
        respx.get(f"{CHEMBL_BASE}/molecule/search").mock(
            return_value=Response(200, json=ASPIRIN_SEARCH_RESPONSE)
        )
        adapter = ChEMBLAdapter()
        results = await adapter.search("aspirin", query_type="name")
        assert len(results) == 1
        compound = results[0]
        assert compound.chembl_id == "CHEMBL25"
        assert compound.name == "ASPIRIN"
        assert "chembl" in compound.sources
        assert compound.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles(self) -> None:
        """Search by SMILES returns Compound objects."""
        respx.get(f"{CHEMBL_BASE}/molecule/search").mock(
            return_value=Response(200, json=ASPIRIN_SEARCH_RESPONSE)
        )
        adapter = ChEMBLAdapter()
        results = await adapter.search("CC(=O)Oc1ccccc1C(O)=O", query_type="smiles")
        assert len(results) >= 1
        assert results[0].chembl_id == "CHEMBL25"

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_chembl_id(self) -> None:
        """Search by ChEMBL ID delegates to get_by_id."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL25").mock(
            return_value=Response(200, json=ASPIRIN_MOLECULE)
        )
        adapter = ChEMBLAdapter()
        results = await adapter.search("CHEMBL25", query_type="identifier")
        assert len(results) == 1
        assert results[0].chembl_id == "CHEMBL25"

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_chembl_prefix_auto_detects_id(self) -> None:
        """Query starting with CHEMBL is auto-detected as ChEMBL ID."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL25").mock(
            return_value=Response(200, json=ASPIRIN_MOLECULE)
        )
        adapter = ChEMBLAdapter()
        # query_type="name" but query starts with "CHEMBL"
        results = await adapter.search("CHEMBL25", query_type="name")
        assert len(results) == 1
        assert results[0].chembl_id == "CHEMBL25"

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_empty_result(self) -> None:
        """Empty molecule list returns empty list."""
        respx.get(f"{CHEMBL_BASE}/molecule/search").mock(
            return_value=Response(200, json={"molecules": [], "page_meta": {}})
        )
        adapter = ChEMBLAdapter()
        results = await adapter.search("nonexistent_compound_xyz")
        assert results == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_not_found_returns_empty_list(self) -> None:
        """search returns empty list when compound is not found."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL99999").mock(
            return_value=Response(404, json={"error": "Not found"})
        )
        adapter = ChEMBLAdapter()
        result = await adapter.search("CHEMBL99999", query_type="identifier")
        assert result == []


class TestChEMBLGetById:
    """Tests for ChEMBL get_by_id."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_by_id_success(self) -> None:
        """get_by_id returns a Compound for a valid ChEMBL ID."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL25").mock(
            return_value=Response(200, json=ASPIRIN_MOLECULE)
        )
        adapter = ChEMBLAdapter()
        compound = await adapter.get_by_id("CHEMBL25")
        assert compound is not None
        assert compound.chembl_id == "CHEMBL25"
        assert compound.name == "ASPIRIN"
        assert compound.smiles == "CC(=O)Oc1ccccc1C(O)=O"
        assert compound.properties.molecular_weight == pytest.approx(180.16)

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_by_id_not_found_returns_none(self) -> None:
        """get_by_id returns None on 404."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL99999").mock(
            return_value=Response(404, json={"error": "Not found"})
        )
        adapter = ChEMBLAdapter()
        compound = await adapter.get_by_id("CHEMBL99999")
        assert compound is None

    @pytest.mark.asyncio
    @respx.mock
    async def test_compound_has_correct_properties(self) -> None:
        """Parsed compound properties are correct."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL25").mock(
            return_value=Response(200, json=ASPIRIN_MOLECULE)
        )
        adapter = ChEMBLAdapter()
        compound = await adapter.get_by_id("CHEMBL25")
        assert compound is not None
        assert compound.properties.xlogp == pytest.approx(1.31)
        assert compound.properties.hbd_count == 1
        assert compound.properties.hba_count == 4
        assert compound.formula == "C9H8O4"


class TestChEMBLBioactivities:
    """Tests for ChEMBL bioactivity retrieval."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_bioactivities_single_page(self) -> None:
        """get_bioactivities returns Bioactivity objects from a single page."""
        single_page = {
            "activities": ASPIRIN_ACTIVITIES_PAGE1["activities"],
            "page_meta": {"next": None, "total_count": 2},
        }
        respx.get(f"{CHEMBL_BASE}/activity").mock(
            return_value=Response(200, json=single_page)
        )
        adapter = ChEMBLAdapter()
        activities = await adapter.get_bioactivities("CHEMBL25")
        assert len(activities) == 2
        assert all(isinstance(a, Bioactivity) for a in activities)
        assert activities[0].activity_type == "IC50"
        assert activities[0].target_name == "Cyclooxygenase-1"
        assert activities[0].value == pytest.approx(1670.0)
        assert activities[0].source == "chembl"

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_bioactivities_pagination(self) -> None:
        """get_bioactivities auto-paginates through multiple pages."""
        # Use side_effect list so the first call returns page 1, second call returns page 2.
        # Both calls hit the same /activity endpoint (with different query params),
        # so a single route with side_effect handles both pages correctly.
        respx.get(f"{CHEMBL_BASE}/activity").mock(
            side_effect=[
                Response(200, json=ASPIRIN_ACTIVITIES_PAGE1),
                Response(200, json=ASPIRIN_ACTIVITIES_PAGE2),
            ]
        )

        adapter = ChEMBLAdapter()
        activities = await adapter.get_bioactivities("CHEMBL25", limit=2)
        assert len(activities) == 3  # 2 from page 1 + 1 from page 2

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_bioactivities_empty(self) -> None:
        """get_bioactivities returns empty list when no activities found."""
        respx.get(f"{CHEMBL_BASE}/activity").mock(
            return_value=Response(200, json={"activities": [], "page_meta": {"next": None}})
        )
        adapter = ChEMBLAdapter()
        activities = await adapter.get_bioactivities("CHEMBL99999")
        assert activities == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_bioactivity_activity_type_normalization(self) -> None:
        """Activity types are normalized by the Bioactivity model."""
        page_data = {
            "activities": [
                {
                    "target_pref_name": "Target A",
                    "standard_type": "ic50",  # lowercase
                    "standard_value": "100",
                    "standard_units": "nM",
                }
            ],
            "page_meta": {"next": None},
        }
        respx.get(f"{CHEMBL_BASE}/activity").mock(
            return_value=Response(200, json=page_data)
        )
        adapter = ChEMBLAdapter()
        activities = await adapter.get_bioactivities("CHEMBL25")
        assert activities[0].activity_type == "IC50"  # normalized

    @pytest.mark.asyncio
    @respx.mock
    async def test_bioactivity_with_none_value(self) -> None:
        """Bioactivity is created even when value is None."""
        page_data = {
            "activities": [
                {
                    "target_pref_name": "Target",
                    "standard_type": "Ki",
                    "standard_value": None,
                }
            ],
            "page_meta": {"next": None},
        }
        respx.get(f"{CHEMBL_BASE}/activity").mock(
            return_value=Response(200, json=page_data)
        )
        adapter = ChEMBLAdapter()
        activities = await adapter.get_bioactivities("CHEMBL25")
        assert len(activities) == 1
        assert activities[0].value is None


class TestChEMBLMechanism:
    """Tests for ChEMBL mechanism of action retrieval."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_mechanism_returns_list(self) -> None:
        """get_mechanism returns a list of mechanism dicts."""
        respx.get(f"{CHEMBL_BASE}/mechanism").mock(
            return_value=Response(200, json=ASPIRIN_MECHANISMS)
        )
        adapter = ChEMBLAdapter()
        mechs = await adapter.get_mechanism("CHEMBL25")
        assert len(mechs) == 1
        assert mechs[0]["action_type"] == "INHIBITOR"
        assert "mechanism_of_action" in mechs[0]

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_mechanism_empty(self) -> None:
        """get_mechanism returns empty list when no mechanisms found."""
        respx.get(f"{CHEMBL_BASE}/mechanism").mock(
            return_value=Response(200, json={"mechanisms": []})
        )
        adapter = ChEMBLAdapter()
        mechs = await adapter.get_mechanism("CHEMBL99999")
        assert mechs == []


class TestChEMBLAdapterMisc:
    """Miscellaneous adapter tests."""

    def test_is_available(self) -> None:
        """is_available always returns True."""
        adapter = ChEMBLAdapter()
        assert adapter.is_available() is True

    def test_adapter_name(self) -> None:
        """Adapter has correct name."""
        adapter = ChEMBLAdapter()
        assert adapter.name == "chembl"

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_properties_returns_dict(self) -> None:
        """get_properties returns a dict."""
        respx.get(f"{CHEMBL_BASE}/molecule/CHEMBL25").mock(
            return_value=Response(200, json=ASPIRIN_MOLECULE)
        )
        adapter = ChEMBLAdapter()
        props = await adapter.get_properties("CHEMBL25")
        assert isinstance(props, dict)
        assert "chembl_id" in props

    @pytest.mark.asyncio
    @respx.mock
    async def test_malformed_response_returns_empty(self) -> None:
        """Malformed API response returns empty list without raising."""
        respx.get(f"{CHEMBL_BASE}/molecule/search").mock(
            return_value=Response(200, json="not a dict")
        )
        adapter = ChEMBLAdapter()
        results = await adapter.search("test")
        assert results == []
