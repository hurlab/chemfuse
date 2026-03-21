"""Integration tests for the complete enrichment pipeline (SPEC-CF-004)."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.models.compound import Compound
from chemfuse.models.patent import Patent
from chemfuse.models.target import TargetAssociation

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
SURECHEMBL_API = "https://www.ebi.ac.uk/surechembl/api/search"
UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/rest"

# ---------------------------------------------------------------------------
# Shared mock payloads
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
ASPIRIN_CHEMBL_ID = "CHEMBL25"

OT_RESPONSE = {
    "data": {
        "drug": {
            "id": ASPIRIN_CHEMBL_ID,
            "name": "ASPIRIN",
            "linkedDiseases": {
                "count": 1,
                "rows": [
                    {
                        "disease": {"id": "EFO_0003785", "name": "pain"},
                        "score": 0.85,
                        "evidenceCount": 42,
                    }
                ],
            },
            "linkedTargets": {
                "count": 1,
                "rows": [
                    {
                        "target": {
                            "id": "ENSG00000095303",
                            "approvedSymbol": "PTGS1",
                            "approvedName": "prostaglandin-endoperoxide synthase 1",
                        }
                    }
                ],
            },
        }
    }
}

SURECHEMBL_RESPONSE = {
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

UNICHEM_RESPONSE = [
    {"src_id": "1", "src_compound_id": ASPIRIN_CHEMBL_ID},
    {"src_id": "22", "src_compound_id": "2244"},
]


class TestFullEnrichmentPipeline:
    """Test complete enrichment pipeline across multiple sources."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_enrich_patents_and_targets(self) -> None:
        """enrich(patents=True, targets=True) populates both fields."""
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))
        respx.get(SURECHEMBL_API).mock(return_value=Response(200, json=SURECHEMBL_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            chembl_id=ASPIRIN_CHEMBL_ID,
            name="aspirin",
        )
        await compound.enrich(patents=True, targets=True)

        assert len(compound.patents) >= 1
        assert len(compound.target_associations) >= 1
        assert all(isinstance(p, Patent) for p in compound.patents)
        assert all(isinstance(a, TargetAssociation) for a in compound.target_associations)

    @pytest.mark.asyncio
    @respx.mock
    async def test_enrich_targets_only(self) -> None:
        """enrich(targets=True) populates target_associations only."""
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            chembl_id=ASPIRIN_CHEMBL_ID,
            name="aspirin",
        )
        await compound.enrich(targets=True)

        assert len(compound.target_associations) >= 1
        assert compound.patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_enrich_patents_only(self) -> None:
        """enrich(patents=True) populates patents only."""
        respx.get(SURECHEMBL_API).mock(return_value=Response(200, json=SURECHEMBL_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            name="aspirin",
        )
        await compound.enrich(patents=True)

        assert len(compound.patents) >= 1
        assert compound.target_associations == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_sources_updated_after_enrichment(self) -> None:
        """compound.sources is updated to include contributing databases."""
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))
        respx.get(SURECHEMBL_API).mock(return_value=Response(200, json=SURECHEMBL_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            chembl_id=ASPIRIN_CHEMBL_ID,
            sources=["pubchem"],
        )
        await compound.enrich(patents=True, targets=True)

        assert "opentargets" in compound.sources
        assert "surechembl" in compound.sources


class TestEnrichmentErrorIsolation:
    """Test that one source failure doesn't block other sources."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_surechembl_down_targets_still_succeed(self) -> None:
        """When SureChEMBL fails, Open Targets results are still returned."""
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))
        # SureChEMBL returns 503
        respx.get(SURECHEMBL_API).mock(return_value=Response(503, text="Service Unavailable"))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            chembl_id=ASPIRIN_CHEMBL_ID,
        )
        # Should not raise an exception
        await compound.enrich(patents=True, targets=True)

        # Targets should be populated
        assert len(compound.target_associations) >= 1
        # Patents should be empty due to failure
        assert compound.patents == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_opentargets_down_patents_still_succeed(self) -> None:
        """When Open Targets fails, SureChEMBL results are still returned."""
        respx.post(OPENTARGETS_URL).mock(return_value=Response(500, text="Error"))
        respx.get(SURECHEMBL_API).mock(return_value=Response(200, json=SURECHEMBL_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            chembl_id=ASPIRIN_CHEMBL_ID,
        )
        await compound.enrich(patents=True, targets=True)

        # Patents should be populated
        assert len(compound.patents) >= 1
        # Targets should be empty due to failure
        assert compound.target_associations == []

    @pytest.mark.asyncio
    async def test_no_exception_on_both_failures(self) -> None:
        """When both sources fail, no exception is raised."""
        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            chembl_id=ASPIRIN_CHEMBL_ID,
        )
        # Both sources will fail with connection errors
        with respx.mock:
            respx.post(OPENTARGETS_URL).mock(return_value=Response(500, text="Error"))
            respx.get(SURECHEMBL_API).mock(return_value=Response(503, text="Error"))
            # Should not raise
            await compound.enrich(patents=True, targets=True)

        assert compound.patents == []
        assert compound.target_associations == []


class TestChemblIdResolution:
    """Test ChEMBL ID resolution via UniChem for compounds without chembl_id."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_resolve_chembl_id_via_unichem_inchikey(self) -> None:
        """Compound with InChIKey resolves ChEMBL ID before querying Open Targets."""
        # UniChem cross-reference by InChIKey
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=UNICHEM_RESPONSE)
        )
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            # No chembl_id initially
        )
        await compound.enrich(targets=True)

        # ChEMBL ID should be resolved and stored
        assert compound.chembl_id == ASPIRIN_CHEMBL_ID
        assert len(compound.target_associations) >= 1

    @pytest.mark.asyncio
    @respx.mock
    async def test_resolve_chembl_id_via_unichem_cid(self) -> None:
        """Compound with CID resolves ChEMBL ID via UniChem when no InChIKey."""
        # UniChem map by PubChem CID
        respx.get(f"{UNICHEM_BASE}/src_compound_id/2244/src_id/22").mock(
            return_value=Response(200, json=UNICHEM_RESPONSE)
        )
        respx.post(OPENTARGETS_URL).mock(return_value=Response(200, json=OT_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            # No inchikey, no chembl_id
        )
        await compound.enrich(targets=True)

        assert compound.chembl_id == ASPIRIN_CHEMBL_ID
        assert len(compound.target_associations) >= 1

    @pytest.mark.asyncio
    @respx.mock
    async def test_skips_opentargets_when_no_chembl_id_resolvable(self) -> None:
        """When ChEMBL ID cannot be resolved, Open Targets is skipped gracefully."""
        # UniChem returns empty
        respx.get(f"{UNICHEM_BASE}/inchikey/{ASPIRIN_INCHIKEY}").mock(
            return_value=Response(200, json=[])
        )
        respx.get(f"{UNICHEM_BASE}/src_compound_id/2244/src_id/22").mock(
            return_value=Response(200, json=[])
        )

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
        )
        # Should not raise even though ChEMBL ID cannot be resolved
        await compound.enrich(targets=True)
        assert compound.target_associations == []


class TestEnrichmentCaching:
    """Test that enrichment respects the skip-if-already-fetched logic."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_skips_source_if_already_in_sources(self) -> None:
        """enrich() skips source if already present in compound.sources."""
        call_count = 0

        def mock_response(request):
            nonlocal call_count
            call_count += 1
            return Response(200, json=OT_RESPONSE)

        respx.post(OPENTARGETS_URL).mock(side_effect=mock_response)
        respx.get(SURECHEMBL_API).mock(return_value=Response(200, json=SURECHEMBL_RESPONSE))

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            inchikey=ASPIRIN_INCHIKEY,
            chembl_id=ASPIRIN_CHEMBL_ID,
            sources=["opentargets"],  # already enriched
        )
        await compound.enrich(targets=True)
        # Should not call Open Targets again (source already present)
        assert call_count == 0

    @pytest.mark.asyncio
    @respx.mock
    async def test_force_refetch_even_if_already_fetched(self) -> None:
        """enrich(force=True) re-fetches even if source is already in sources."""
        call_count = 0

        def mock_response(request):
            nonlocal call_count
            call_count += 1
            return Response(200, json=OT_RESPONSE)

        respx.post(OPENTARGETS_URL).mock(side_effect=mock_response)

        compound = Compound(
            cid=2244,
            smiles=ASPIRIN_SMILES,
            chembl_id=ASPIRIN_CHEMBL_ID,
            sources=["opentargets"],  # already enriched
        )
        await compound.enrich(targets=True, force=True)
        # With force=True, should still call Open Targets
        assert call_count == 1


class TestEnrichmentExtendedField:
    """Tests for the patents and target_associations fields on Compound."""

    def test_compound_has_target_associations_field(self) -> None:
        """Compound has target_associations field defaulting to empty list."""
        compound = Compound(smiles="CC")
        assert compound.target_associations == []

    def test_compound_has_patents_field(self) -> None:
        """Compound has patents field defaulting to empty list."""
        compound = Compound(smiles="CC")
        assert compound.patents == []

    def test_compound_to_dict_includes_new_fields(self) -> None:
        """Compound.to_dict() includes target_associations and patents."""
        compound = Compound(smiles="CC")
        d = compound.to_dict()
        assert "target_associations" in d
        assert "patents" in d
