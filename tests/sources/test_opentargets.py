"""Tests for Open Targets GraphQL source adapter."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.core.exceptions import SourceError
from chemfuse.models.target import TargetAssociation
from chemfuse.sources.opentargets import OpenTargetsAdapter

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# ---------------------------------------------------------------------------
# Mock response fixtures
# ---------------------------------------------------------------------------

ASPIRIN_CHEMBL_ID = "CHEMBL25"

OPENTARGETS_DISEASE_RESPONSE = {
    "data": {
        "drug": {
            "id": ASPIRIN_CHEMBL_ID,
            "name": "ASPIRIN",
            "linkedDiseases": {
                "count": 1,
                "rows": [
                    {
                        "disease": {
                            "id": "EFO_0003785",
                            "name": "pain",
                        },
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

OPENTARGETS_EMPTY_RESPONSE = {
    "data": {
        "drug": None,
    }
}

OPENTARGETS_ERROR_RESPONSE = {
    "errors": [
        {
            "message": "Variable 'chemblId' not found",
            "locations": [{"line": 2, "column": 9}],
        }
    ]
}

OPENTARGETS_MULTIPLE_DISEASES_RESPONSE = {
    "data": {
        "drug": {
            "id": "CHEMBL25",
            "name": "ASPIRIN",
            "linkedDiseases": {
                "count": 2,
                "rows": [
                    {
                        "disease": {"id": "EFO_001", "name": "inflammation"},
                        "score": 0.9,
                        "evidenceCount": 100,
                    },
                    {
                        "disease": {"id": "EFO_002", "name": "fever"},
                        "score": 0.6,
                        "evidenceCount": 50,
                    },
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

# Response with 5 diseases and 3 targets — used to verify no Cartesian product
OPENTARGETS_5D_3T_RESPONSE = {
    "data": {
        "drug": {
            "id": "CHEMBL42",
            "name": "TESTDRUG",
            "linkedDiseases": {
                "count": 5,
                "rows": [
                    {"disease": {"id": f"EFO_00{i}", "name": f"disease_{i}"}, "score": 0.5, "evidenceCount": 10}
                    for i in range(1, 6)
                ],
            },
            "linkedTargets": {
                "count": 3,
                "rows": [
                    {"target": {"id": f"ENSG0000000000{i}", "approvedSymbol": f"GENE{i}", "approvedName": f"gene {i}"}}
                    for i in range(1, 4)
                ],
            },
        }
    }
}

OPENTARGETS_NO_DISEASES_RESPONSE = {
    "data": {
        "drug": {
            "id": "CHEMBL999",
            "name": "NOVEL",
            "linkedDiseases": {
                "count": 0,
                "rows": [],
            },
            "linkedTargets": {
                "count": 0,
                "rows": [],
            },
        }
    }
}

OPENTARGETS_TARGET_INFO_RESPONSE = {
    "data": {
        "target": {
            "id": "ENSG00000095303",
            "approvedSymbol": "PTGS1",
            "tractability": [
                {"label": "High Quality Ligand", "modality": "SM", "value": True}
            ],
        }
    }
}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestOpenTargetsSearchByChemblId:
    """Tests for OpenTargetsAdapter.search_by_chembl_id."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_returns_target_associations(self) -> None:
        """search_by_chembl_id returns TargetAssociation objects for aspirin."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_DISEASE_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)
        assert len(associations) >= 1
        assert all(isinstance(a, TargetAssociation) for a in associations)

    @pytest.mark.asyncio
    @respx.mock
    async def test_association_fields_populated(self) -> None:
        """search_by_chembl_id correctly populates TargetAssociation fields.

        With one disease and one target the response produces 2 independent
        associations: one disease association and one target association.
        """
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_DISEASE_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)
        # Exactly 2 associations: one disease, one target (no Cartesian product)
        assert len(associations) == 2

        disease_assoc = next(a for a in associations if a.disease_name == "pain")
        assert disease_assoc.association_score == pytest.approx(0.85)
        assert disease_assoc.evidence_count == 42

        target_assoc = next(a for a in associations if a.target_name == "PTGS1")
        assert target_assoc.disease_name is None

    @pytest.mark.asyncio
    @respx.mock
    async def test_source_is_opentargets(self) -> None:
        """All returned associations have source='opentargets'."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_DISEASE_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)
        assert all(a.source == "opentargets" for a in associations)

    @pytest.mark.asyncio
    @respx.mock
    async def test_empty_result_when_drug_not_found(self) -> None:
        """search_by_chembl_id returns empty list when drug is not in Open Targets."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_EMPTY_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id("CHEMBL99999999")
        assert associations == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_empty_result_when_no_diseases(self) -> None:
        """search_by_chembl_id returns empty list when drug has no disease associations."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_NO_DISEASES_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id("CHEMBL999")
        assert associations == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_graphql_errors_raise_source_error(self) -> None:
        """search_by_chembl_id raises SourceError when GraphQL response contains errors."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_ERROR_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        with pytest.raises(SourceError) as exc_info:
            await adapter.search_by_chembl_id("BADID")
        assert "GraphQL" in str(exc_info.value) or "Variable" in str(exc_info.value)

    @pytest.mark.asyncio
    @respx.mock
    async def test_multiple_diseases_returns_multiple_associations(self) -> None:
        """search_by_chembl_id returns one association per disease plus one per target.

        2 diseases + 1 target must produce exactly 3 associations, not 2*1=2 or
        a Cartesian 2*1=2. The two types are emitted independently.
        """
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_MULTIPLE_DISEASES_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)
        # 2 disease associations + 1 target association = 3 total
        assert len(associations) == 3
        disease_names = [a.disease_name for a in associations if a.disease_name]
        assert "inflammation" in disease_names
        assert "fever" in disease_names
        target_names = [a.target_name for a in associations]
        assert "PTGS1" in target_names

    @pytest.mark.asyncio
    @respx.mock
    async def test_association_score_normalized(self) -> None:
        """Association scores are normalized to [0, 1] range."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_DISEASE_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)
        for assoc in associations:
            if assoc.association_score is not None:
                assert 0.0 <= assoc.association_score <= 1.0

    @pytest.mark.asyncio
    @respx.mock
    async def test_api_unavailable_raises_source_error(self) -> None:
        """HTTP 500 from Open Targets API raises an exception that propagates."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(500, text="Internal Server Error")
        )
        adapter = OpenTargetsAdapter()
        # Post method raises SourceError on non-200
        with pytest.raises(SourceError):
            await adapter.search_by_chembl_id(ASPIRIN_CHEMBL_ID)

    @pytest.mark.asyncio
    @respx.mock
    async def test_5_diseases_3_targets_produces_8_not_15(self) -> None:
        """5 diseases + 3 targets must produce exactly 8 associations, not 15.

        The disease and target lists are independent; they must not be
        cross-joined into a Cartesian product (5*3=15).
        """
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_5D_3T_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id("CHEMBL42")
        assert len(associations) == 8, (
            f"Expected 8 associations (5 disease + 3 target), got {len(associations)}. "
            "Disease and target associations must not be cross-joined."
        )

    @pytest.mark.asyncio
    @respx.mock
    async def test_disease_associations_have_no_target_cross_product(self) -> None:
        """Disease associations carry disease info; target fields are not replicated per disease."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_5D_3T_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id("CHEMBL42")
        disease_assocs = [a for a in associations if a.disease_name is not None]
        # Exactly 5 disease associations, one per disease row
        assert len(disease_assocs) == 5
        disease_names = {a.disease_name for a in disease_assocs}
        assert disease_names == {f"disease_{i}" for i in range(1, 6)}

    @pytest.mark.asyncio
    @respx.mock
    async def test_target_associations_are_independent(self) -> None:
        """Target associations carry target info; disease fields are None."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_5D_3T_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        associations = await adapter.search_by_chembl_id("CHEMBL42")
        target_assocs = [a for a in associations if a.disease_name is None]
        # Exactly 3 target associations, one per target row
        assert len(target_assocs) == 3
        target_names = {a.target_name for a in target_assocs}
        assert target_names == {"GENE1", "GENE2", "GENE3"}
        # Confirm disease fields are absent on target-only entries
        for ta in target_assocs:
            assert ta.disease_id is None
            assert ta.disease_name is None


class TestOpenTargetsGetTargetInfo:
    """Tests for OpenTargetsAdapter.get_target_info."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_target_info_returns_dict(self) -> None:
        """get_target_info returns a dict with target information."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_TARGET_INFO_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        info = await adapter.get_target_info("ENSG00000095303")
        assert isinstance(info, dict)
        assert info.get("id") == "ENSG00000095303"
        assert info.get("approvedSymbol") == "PTGS1"

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_target_info_graphql_error_raises(self) -> None:
        """get_target_info raises SourceError on GraphQL errors."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json=OPENTARGETS_ERROR_RESPONSE)
        )
        adapter = OpenTargetsAdapter()
        with pytest.raises(SourceError):
            await adapter.get_target_info("ENSG_INVALID")

    @pytest.mark.asyncio
    @respx.mock
    async def test_get_target_info_empty_returns_empty_dict(self) -> None:
        """get_target_info returns empty dict when target not found."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json={"data": {"target": None}})
        )
        adapter = OpenTargetsAdapter()
        info = await adapter.get_target_info("ENSG_NOT_FOUND")
        assert info == {}


class TestOpenTargetsAdapterMisc:
    """Miscellaneous OpenTargets adapter tests."""

    def test_adapter_name(self) -> None:
        """Adapter has correct name."""
        adapter = OpenTargetsAdapter()
        assert adapter.name == "opentargets"

    def test_is_available(self) -> None:
        """is_available always returns True."""
        adapter = OpenTargetsAdapter()
        assert adapter.is_available() is True

    @pytest.mark.asyncio
    async def test_search_returns_empty(self) -> None:
        """search always returns empty list."""
        adapter = OpenTargetsAdapter()
        result = await adapter.search("aspirin")
        assert result == []

    @pytest.mark.asyncio
    async def test_get_by_id_returns_none(self) -> None:
        """get_by_id always returns None."""
        adapter = OpenTargetsAdapter()
        result = await adapter.get_by_id("CHEMBL25")
        assert result is None

    @pytest.mark.asyncio
    async def test_get_properties_returns_empty(self) -> None:
        """get_properties always returns empty dict."""
        adapter = OpenTargetsAdapter()
        result = await adapter.get_properties("CHEMBL25")
        assert result == {}

    def test_rate_limit_attribute(self) -> None:
        """Adapter has rate_limit attribute."""
        adapter = OpenTargetsAdapter()
        assert adapter.rate_limit > 0

    @pytest.mark.asyncio
    @respx.mock
    async def test_parse_response_with_no_data_key(self) -> None:
        """Response with no 'data' key returns empty list."""
        respx.post(OPENTARGETS_URL).mock(
            return_value=Response(200, json={})
        )
        adapter = OpenTargetsAdapter()
        # _parse_drug_disease_response should return [] for empty response
        result = adapter._parse_drug_disease_response({}, "CHEMBL25")
        assert result == []

    def test_parse_response_with_errors(self) -> None:
        """_parse_drug_disease_response raises SourceError on errors field."""
        adapter = OpenTargetsAdapter()
        with pytest.raises(SourceError):
            adapter._parse_drug_disease_response(OPENTARGETS_ERROR_RESPONSE, "CHEMBL25")

    def test_parse_response_non_dict_returns_empty(self) -> None:
        """_parse_drug_disease_response returns [] for non-dict input."""
        adapter = OpenTargetsAdapter()
        result = adapter._parse_drug_disease_response("not a dict", "CHEMBL25")
        assert result == []
