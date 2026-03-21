"""Integration tests for multi-source search, enrichment, and cross-reference."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest
import respx
from httpx import Response

import chemfuse as cf
from chemfuse.models.bioactivity import BindingMeasurement, Bioactivity
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/rest"
BINDINGDB_BASE = "https://www.bindingdb.org/axis2/services/BDBService"

ASPIRIN_INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

# ChEMBL aspirin molecule
CHEMBL_ASPIRIN = {
    "molecule_chembl_id": "CHEMBL25",
    "pref_name": "ASPIRIN",
    "molecule_structures": {
        "canonical_smiles": "CC(=O)Oc1ccccc1C(O)=O",
        "standard_inchi_key": ASPIRIN_INCHIKEY,
    },
    "molecule_properties": {
        "full_mwt": "180.16",
        "alogp": "1.31",
        "molecular_formula": "C9H8O4",
    },
}

CHEMBL_SEARCH_RESPONSE = {
    "molecules": [CHEMBL_ASPIRIN],
    "page_meta": {"next": None, "total_count": 1},
}

CHEMBL_ACTIVITIES_RESPONSE = {
    "activities": [
        {
            "target_pref_name": "Cyclooxygenase-1",
            "target_chembl_id": "CHEMBL247",
            "standard_type": "IC50",
            "standard_value": "1670",
            "standard_units": "nM",
            "standard_relation": "=",
        }
    ],
    "page_meta": {"next": None},
}

UNICHEM_ASPIRIN = [
    {"src_id": "1", "src_compound_id": "CHEMBL25"},
    {"src_id": "22", "src_compound_id": "2244"},
    {"src_id": "2", "src_compound_id": "DB00945"},
]

BINDINGDB_ASPIRIN = {
    "affinities": [
        {"target": "Cyclooxygenase-1", "uniprot": "P23219", "ki": "1670"},
    ]
}


# --- Helpers ---

def make_pubchem_aspirin() -> Compound:
    """Create a Compound representing aspirin from PubChem."""
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey=ASPIRIN_INCHIKEY,
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem"],
        properties=CompoundProperties(molecular_weight=180.16),
    )


def make_chembl_aspirin() -> Compound:
    """Create a Compound representing aspirin from ChEMBL."""
    return Compound(
        chembl_id="CHEMBL25",
        smiles="CC(=O)Oc1ccccc1C(O)=O",
        inchikey=ASPIRIN_INCHIKEY,
        name="ASPIRIN",
        sources=["chembl"],
        properties=CompoundProperties(molecular_weight=180.16, xlogp=1.31),
    )


# --- Multi-source search tests ---

class TestMultiSourceSearch:
    """Tests for parallel multi-source search and InChIKey-based merge."""

    @pytest.mark.asyncio
    async def test_search_merges_by_inchikey(self) -> None:
        """Multi-source search merges compounds with same InChIKey."""
        pubchem_compound = make_pubchem_aspirin()
        chembl_compound = make_chembl_aspirin()

        with (
            patch("chemfuse.registry") as mock_registry,
        ):
            mock_pubchem = AsyncMock()
            mock_pubchem.search = AsyncMock(return_value=[pubchem_compound])
            mock_chembl = AsyncMock()
            mock_chembl.search = AsyncMock(return_value=[chembl_compound])

            def get_adapter(name: str):
                if name == "pubchem":
                    return mock_pubchem
                return mock_chembl

            mock_registry.get.side_effect = get_adapter

            collection = await cf.search_async("aspirin", sources=["pubchem", "chembl"])

        # Should have exactly one compound (merged)
        assert len(collection) == 1
        merged = collection[0]
        # Has identifiers from both sources
        assert merged.cid == 2244
        assert merged.chembl_id == "CHEMBL25"
        # Sources list includes both
        assert "pubchem" in merged.sources
        assert "chembl" in merged.sources

    @pytest.mark.asyncio
    async def test_search_source_failure_returns_partial_results(self) -> None:
        """When one source fails, results from other sources are returned."""
        pubchem_compound = make_pubchem_aspirin()

        with patch("chemfuse.registry") as mock_registry:
            mock_pubchem = AsyncMock()
            mock_pubchem.search = AsyncMock(return_value=[pubchem_compound])

            mock_chembl = AsyncMock()
            mock_chembl.search = AsyncMock(side_effect=Exception("ChEMBL is down"))

            def get_adapter(name: str):
                if name == "pubchem":
                    return mock_pubchem
                return mock_chembl

            mock_registry.get.side_effect = get_adapter

            collection = await cf.search_async("aspirin", sources=["pubchem", "chembl"])

        # PubChem results still returned
        assert len(collection) >= 1
        assert "pubchem" in collection.sources
        # ChEMBL not in successful sources
        assert "chembl" not in collection.sources

    @pytest.mark.asyncio
    async def test_search_progress_callback_called(self) -> None:
        """Progress callback is called for each source."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            called_sources = []

            def progress(source_name: str) -> None:
                called_sources.append(source_name)

            await cf.search_async(
                "aspirin",
                sources=["pubchem", "chembl"],
                progress_callback=progress,
            )

        assert "pubchem" in called_sources
        assert "chembl" in called_sources

    @pytest.mark.asyncio
    async def test_search_no_inchikey_compounds_not_merged(self) -> None:
        """Compounds without InChIKey are not merged."""
        c1 = Compound(smiles="CC", sources=["pubchem"])
        c2 = Compound(smiles="CC", sources=["chembl"])

        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(side_effect=[
                [c1],  # pubchem call
                [c2],  # chembl call
            ])
            mock_registry.get.return_value = mock_adapter

            collection = await cf.search_async("test", sources=["pubchem", "chembl"])

        # Both kept since no InChIKey to deduplicate on
        assert len(collection) == 2


# --- Cross-reference tests ---

class TestCrossReference:
    """Tests for cf.map_identifiers cross-reference function."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_identifiers_by_cid(self) -> None:
        """map_identifiers_async maps PubChem CID to other databases."""
        respx.get(f"{UNICHEM_BASE}/src_compound_id/2244/src_id/22").mock(
            return_value=Response(200, json=UNICHEM_ASPIRIN)
        )
        result = await cf.map_identifiers_async(cid=2244)
        assert result["chembl"] == "CHEMBL25"
        assert result["pubchem"] == "2244"
        assert result["drugbank"] == "DB00945"

    @pytest.mark.asyncio
    @respx.mock
    async def test_map_identifiers_by_inchikey(self) -> None:
        """map_identifiers_async maps InChIKey to all databases."""
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=UNICHEM_ASPIRIN)
        )
        result = await cf.map_identifiers_async(inchikey=ASPIRIN_INCHIKEY)
        assert "chembl" in result
        assert "pubchem" in result

    @pytest.mark.asyncio
    async def test_map_identifiers_returns_input_on_empty_result(self) -> None:
        """map_identifiers returns at least the input if UniChem finds nothing."""
        with patch("chemfuse.sources.unichem.UniChemAdapter.map_identifiers", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {}
            result = await cf.map_identifiers_async(cid=2244)

        assert result.get("pubchem") == "2244"

    @pytest.mark.asyncio
    async def test_map_identifiers_no_args_returns_empty_like_result(self) -> None:
        """map_identifiers with no useful args returns a dict."""
        result = await cf.map_identifiers_async()
        # Returns empty or minimal dict
        assert isinstance(result, dict)


# --- Enrichment tests ---

class TestCompoundEnrichment:
    """Tests for Compound.enrich() method."""

    @pytest.mark.asyncio
    async def test_enrich_chembl_populates_bioactivities(self) -> None:
        """enrich with chembl source fetches and stores bioactivities."""
        compound = Compound(
            cid=2244,
            chembl_id="CHEMBL25",
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey=ASPIRIN_INCHIKEY,
            name="aspirin",
            sources=["pubchem"],
        )

        mock_activities = [
            Bioactivity(
                target_name="COX-1",
                activity_type="IC50",
                value=1670.0,
                units="nM",
                source="chembl",
            )
        ]

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_chembl_adapter = AsyncMock()
            mock_chembl_adapter.get_bioactivities = AsyncMock(return_value=mock_activities)
            mock_registry.get.return_value = mock_chembl_adapter

            await compound.enrich(sources=["chembl"])

        assert len(compound.bioactivities) == 1
        assert compound.bioactivities[0].target_name == "COX-1"
        assert "chembl" in compound.sources

    @pytest.mark.asyncio
    async def test_enrich_binding_populates_binding_data(self) -> None:
        """enrich with binding=True fetches binding measurements."""
        compound = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey=ASPIRIN_INCHIKEY,
            name="aspirin",
            sources=["pubchem"],
        )

        mock_measurements = [
            BindingMeasurement(
                target_name="COX-1",
                target_uniprot="P23219",
                ki=1670.0,
            )
        ]

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_bindingdb = AsyncMock()
            mock_bindingdb.search_by_smiles = AsyncMock(return_value=mock_measurements)
            mock_registry.get.return_value = mock_bindingdb

            await compound.enrich(binding=True)

        assert len(compound.binding_data) == 1
        assert compound.binding_data[0].ki == pytest.approx(1670.0)
        assert "bindingdb" in compound.sources

    @pytest.mark.asyncio
    async def test_enrich_lazy_skips_existing_sources(self) -> None:
        """enrich skips sources already in compound.sources."""
        compound = Compound(
            cid=2244,
            smiles="CC",
            sources=["pubchem", "chembl"],  # chembl already fetched
        )

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_bioactivities = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            await compound.enrich(sources=["chembl"])

        # get_bioactivities should NOT have been called (lazy loading)
        mock_registry.get.assert_not_called()

    @pytest.mark.asyncio
    async def test_enrich_force_refetches(self) -> None:
        """enrich with force=True re-fetches even if source already present."""
        compound = Compound(
            cid=2244,
            chembl_id="CHEMBL25",
            smiles="CC",
            sources=["pubchem", "chembl"],
        )

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_bioactivities = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            await compound.enrich(sources=["chembl"], force=True)

        # Should have been called despite chembl already in sources
        mock_registry.get.assert_called()

    @pytest.mark.asyncio
    async def test_enrich_failure_logs_warning_not_raises(self) -> None:
        """enrich failure logs a warning but does not raise."""
        compound = Compound(
            cid=2244,
            chembl_id="CHEMBL25",
            smiles="CC",
            sources=["pubchem"],
        )

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_bioactivities = AsyncMock(
                side_effect=Exception("Network error")
            )
            mock_registry.get.return_value = mock_adapter

            # Should not raise
            await compound.enrich(sources=["chembl"])

        assert compound.bioactivities == []

    @pytest.mark.asyncio
    async def test_enrich_does_not_duplicate_binding_on_repeat(self) -> None:
        """Calling enrich(binding=True) twice with force=False skips second fetch."""
        compound = Compound(
            cid=2244,
            smiles="CC",
            sources=["pubchem", "bindingdb"],  # already enriched
        )

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search_by_smiles = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            await compound.enrich(binding=True)

        # Should NOT call since bindingdb already in sources
        mock_registry.get.assert_not_called()


# --- CompoundCollection.enrich_all tests ---

class TestCollectionEnrichAll:
    """Tests for CompoundCollection.enrich_all() method."""

    @pytest.mark.asyncio
    async def test_enrich_all_calls_enrich_on_each_compound(self) -> None:
        """enrich_all calls enrich() on every compound."""
        from chemfuse.models.compound import Compound as CompoundClass

        compounds = [
            Compound(smiles="CC", name="methane", sources=["pubchem"]),
            Compound(smiles="CCC", name="propane", sources=["pubchem"]),
        ]
        collection = CompoundCollection(compounds=compounds, sources=["pubchem"])

        enrich_calls = []

        async def mock_enrich(self, sources=None, binding=False, force=False):
            enrich_calls.append(1)

        with patch.object(CompoundClass, "enrich", mock_enrich):
            await collection.enrich_all(sources=["chembl"])

        assert len(enrich_calls) == 2

    @pytest.mark.asyncio
    async def test_enrich_all_progress_callback(self) -> None:
        """enrich_all calls progress_callback with correct arguments."""
        from chemfuse.models.compound import Compound as CompoundClass

        compound = Compound(smiles="CC", name="methane", sources=["pubchem"])
        collection = CompoundCollection(compounds=[compound], sources=["pubchem"])

        progress_records = []

        def progress(completed: int, total: int, compound_name: str) -> None:
            progress_records.append((completed, total, compound_name))

        async def mock_enrich(self, sources=None, binding=False, force=False):
            pass

        with patch.object(CompoundClass, "enrich", mock_enrich):
            await collection.enrich_all(sources=["chembl"], progress_callback=progress)

        assert len(progress_records) == 1
        completed, total, name = progress_records[0]
        assert completed == 1
        assert total == 1
        assert name == "methane"

    @pytest.mark.asyncio
    async def test_enrich_all_continues_on_failure(self) -> None:
        """enrich_all continues when individual compound enrichment fails."""
        from chemfuse.models.compound import Compound as CompoundClass

        c1 = Compound(smiles="CC", name="methane", sources=["pubchem"])
        c2 = Compound(smiles="CCC", name="propane", sources=["pubchem"])
        collection = CompoundCollection(compounds=[c1, c2], sources=["pubchem"])

        enrich_count = [0]
        call_idx = [0]

        async def mock_enrich_dispatch(self, sources=None, binding=False, force=False):
            idx = call_idx[0]
            call_idx[0] += 1
            if idx == 0:
                raise Exception("Failed")
            enrich_count[0] += 1

        with patch.object(CompoundClass, "enrich", mock_enrich_dispatch):
            # Should not raise
            await collection.enrich_all(sources=["chembl"])

        assert enrich_count[0] == 1  # c2 succeeded
