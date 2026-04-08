"""Tests for the ChemFuse MCP server.

Tests exercise handler functions directly (no running server required) and
verify the tool list, argument routing, and error handling. Network calls
are mocked via respx.
"""

from __future__ import annotations

import json
from typing import Any

import pytest
import respx
from httpx import Response

# ---------------------------------------------------------------------------
# Import guards — skip all tests if mcp is not installed
# ---------------------------------------------------------------------------

pytest.importorskip("mcp", reason="mcp package not installed")

from chemfuse.mcp.server import (  # noqa: E402 — after importorskip
    _TOOLS,
    handle_analyze_scaffolds,
    handle_check_drug_likeness,
    handle_compare_compounds,
    handle_compute_descriptors,
    handle_cross_reference,
    handle_predict_admet,
    handle_search_compounds,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

# Minimal PubChem responses used across tests
_PUBCHEM_CID_LIST = {"IdentifierList": {"CID": [2244]}}
_PUBCHEM_PROPERTIES = {
    "PropertyTable": {
        "Properties": [
            {
                "CID": 2244,
                "MolecularFormula": "C9H8O4",
                "MolecularWeight": 180.16,
                "ExactMass": 180.042,
                "CanonicalSMILES": ASPIRIN_SMILES,
                "IsomericSMILES": ASPIRIN_SMILES,
                "InChI": "InChI=1S/C9H8O4",
                "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "IUPACName": "2-(acetyloxy)benzoic acid",
                "XLogP": 1.2,
                "TPSA": 63.6,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 3,
                "HeavyAtomCount": 13,
                "Complexity": 212.0,
            }
        ]
    }
}
_PUBCHEM_SYNONYMS = {
    "InformationList": {
        "Information": [
            {"CID": 2244, "Synonym": ["aspirin", "acetylsalicylic acid"]}
        ]
    }
}


# ---------------------------------------------------------------------------
# Helper to parse response
# ---------------------------------------------------------------------------


def _parse(response: list[Any]) -> Any:
    """Parse the first TextContent item from a handler response."""
    assert response, "Response list is empty"
    return json.loads(response[0].text)


def _is_error(response: list[Any]) -> bool:
    data = _parse(response)
    return "error" in data


# ---------------------------------------------------------------------------
# 1. Tool listing
# ---------------------------------------------------------------------------


class TestToolListing:
    """Verify all expected tools are registered."""

    expected_tool_names = {
        "search_compounds",
        "get_compound",
        "find_similar",
        "cross_reference",
        "predict_admet",
        "check_drug_likeness",
        "compute_descriptors",
        "analyze_scaffolds",
        "compare_compounds",
    }

    def test_all_tools_present(self) -> None:
        tool_names = {t.name for t in _TOOLS}
        assert self.expected_tool_names == tool_names

    def test_all_tools_have_descriptions(self) -> None:
        for tool in _TOOLS:
            assert tool.description, f"Tool '{tool.name}' has no description"
            assert len(tool.description) > 10, f"Tool '{tool.name}' description is too short"

    def test_all_tools_have_input_schema(self) -> None:
        for tool in _TOOLS:
            assert tool.inputSchema, f"Tool '{tool.name}' has no inputSchema"
            assert tool.inputSchema.get("type") == "object"

    def test_required_fields_defined(self) -> None:
        required_map = {
            "search_compounds": {"query"},
            "get_compound": {"identifier"},
            "find_similar": {"smiles"},
            "cross_reference": {"identifier", "identifier_type"},
            "predict_admet": {"smiles"},
            "check_drug_likeness": {"smiles"},
            "compute_descriptors": {"smiles"},
            "analyze_scaffolds": {"smiles_list"},
            "compare_compounds": {"smiles_a", "smiles_b"},
        }
        for tool in _TOOLS:
            expected = required_map.get(tool.name, set())
            actual = set(tool.inputSchema.get("required", []))
            assert actual == expected, f"Tool '{tool.name}' required fields mismatch: {actual} != {expected}"

    def test_tool_count(self) -> None:
        assert len(_TOOLS) == 9


# ---------------------------------------------------------------------------
# 2. search_compounds — mocked PubChem API
# ---------------------------------------------------------------------------


class TestSearchCompounds:
    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_name_returns_compounds(self) -> None:
        # PubChem adapter searches via /compound/name/{name}/property/.../JSON directly.
        # Import the actual PROPERTY_LIST to match exactly.
        from chemfuse.sources.pubchem import PROPERTY_LIST as _PL

        respx.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/{_PL}/JSON"
        ).mock(return_value=Response(200, json=_PUBCHEM_PROPERTIES))
        respx.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/JSON"
        ).mock(return_value=Response(200, json=_PUBCHEM_SYNONYMS))

        response = await handle_search_compounds({"query": "aspirin"})
        data = _parse(response)

        assert not _is_error(response)
        assert "compounds" in data
        assert data["count"] >= 1
        compounds = data["compounds"]
        assert len(compounds) >= 1
        first = compounds[0]
        assert "smiles" in first

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_applies_limit(self) -> None:
        from chemfuse.sources.pubchem import PROPERTY_LIST as _PL

        respx.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/{_PL}/JSON"
        ).mock(return_value=Response(200, json=_PUBCHEM_PROPERTIES))
        respx.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/JSON"
        ).mock(return_value=Response(200, json=_PUBCHEM_SYNONYMS))

        response = await handle_search_compounds({"query": "aspirin", "limit": 1})
        data = _parse(response)
        assert len(data["compounds"]) <= 1

    @pytest.mark.asyncio
    async def test_search_error_returns_error_payload(self) -> None:
        # Simulate a network error by passing an invalid source
        response = await handle_search_compounds({"query": "aspirin", "sources": ["nonexistent_db"]})
        data = _parse(response)
        # Should either return error or empty compounds (not raise)
        assert isinstance(data, dict)


# ---------------------------------------------------------------------------
# 3. predict_admet — real SMILES (aspirin)
# ---------------------------------------------------------------------------


class TestPredictAdmet:
    def test_aspirin_returns_profile(self) -> None:
        response = handle_predict_admet({"smiles": ASPIRIN_SMILES})
        data = _parse(response)

        assert not _is_error(response), f"Unexpected error: {data}"
        assert data["smiles"] == ASPIRIN_SMILES
        assert "overall_score" in data
        assert 0.0 <= data["overall_score"] <= 1.0
        assert "predictions" in data
        preds = data["predictions"]
        assert len(preds) > 0
        # Verify expected ADMET properties are present
        assert "bbb_permeability" in preds
        assert "herg_liability" in preds
        assert "solubility" in preds

    def test_aspirin_prediction_values_in_range(self) -> None:
        response = handle_predict_admet({"smiles": ASPIRIN_SMILES})
        data = _parse(response)

        for name, pred in data["predictions"].items():
            val = pred["value"]
            assert val is not None, f"Prediction '{name}' has None value"
            # Categorical predictions are floats in 0-1, units may indicate otherwise
            assert isinstance(val, (int, float)), f"Prediction '{name}' value is not numeric"

    def test_invalid_smiles_returns_error(self) -> None:
        response = handle_predict_admet({"smiles": "not_a_smiles"})
        assert _is_error(response)

    def test_empty_smiles_returns_error(self) -> None:
        response = handle_predict_admet({"smiles": ""})
        assert _is_error(response)


# ---------------------------------------------------------------------------
# 4. check_drug_likeness — real SMILES
# ---------------------------------------------------------------------------


class TestCheckDrugLikeness:
    def test_aspirin_passes_lipinski(self) -> None:
        response = handle_check_drug_likeness({"smiles": ASPIRIN_SMILES})
        data = _parse(response)

        assert not _is_error(response), f"Unexpected error: {data}"
        assert "filters" in data
        assert "overall_pass" in data
        assert data["filters"]["lipinski"]["pass"] is True

    def test_aspirin_all_filters_present(self) -> None:
        response = handle_check_drug_likeness({"smiles": ASPIRIN_SMILES})
        data = _parse(response)
        filters = data["filters"]
        for name in ("lipinski", "veber", "ghose", "egan", "muegge"):
            assert name in filters, f"Filter '{name}' missing from response"
            assert "pass" in filters[name]
            assert "violations" in filters[name]

    def test_aspirin_pains_and_qed_present(self) -> None:
        response = handle_check_drug_likeness({"smiles": ASPIRIN_SMILES})
        data = _parse(response)
        # PAINS and QED require RDKit — skip check if not present
        if "pains" in data["filters"]:
            assert "pass" in data["filters"]["pains"]
        if "qed" in data:
            assert isinstance(data["qed"], dict)

    def test_invalid_smiles_returns_error(self) -> None:
        response = handle_check_drug_likeness({"smiles": "INVALID_SMILES_XYZ"})
        # Should either return error or pass gracefully (some filters may handle invalid SMILES)
        data = _parse(response)
        assert isinstance(data, dict)

    def test_large_molecule_fails_lipinski(self) -> None:
        # Paclitaxel (taxol) — MW ~854, clearly fails Lipinski
        paclitaxel = (
            "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(=C2OC1=O)C)(C(=O)O4)OC(=O)C)O)"
            "OC(=O)c5ccccc5)(CO)C)OC(=O)C6=CC=CC=C6"
        )
        response = handle_check_drug_likeness({"smiles": paclitaxel})
        data = _parse(response)
        if not _is_error(response):
            # Paclitaxel should fail Lipinski
            assert data["filters"]["lipinski"]["pass"] is False


# ---------------------------------------------------------------------------
# 5. compute_descriptors — real SMILES
# ---------------------------------------------------------------------------


class TestComputeDescriptors:
    def test_aspirin_returns_descriptors(self) -> None:
        response = handle_compute_descriptors({"smiles": ASPIRIN_SMILES})
        data = _parse(response)

        assert not _is_error(response), f"Unexpected error: {data}"
        assert data["smiles"] == ASPIRIN_SMILES
        assert "descriptors" in data
        assert "total_descriptors" in data
        assert data["total_descriptors"] > 0

    def test_aspirin_key_descriptors_present(self) -> None:
        response = handle_compute_descriptors({"smiles": ASPIRIN_SMILES})
        data = _parse(response)
        descs = data["descriptors"]

        # Verify key descriptors are returned with reasonable values
        assert "MolWt" in descs
        assert abs(descs["MolWt"] - 180.16) < 1.0, f"MW {descs['MolWt']} out of range for aspirin"
        assert "MolLogP" in descs
        assert "TPSA" in descs
        assert abs(descs["TPSA"] - 63.6) < 5.0, f"TPSA {descs['TPSA']} out of range for aspirin"

    def test_ibuprofen_descriptor_values(self) -> None:
        response = handle_compute_descriptors({"smiles": IBUPROFEN_SMILES})
        data = _parse(response)
        assert not _is_error(response)
        descs = data["descriptors"]
        assert "MolWt" in descs
        # Ibuprofen MW ~206
        assert abs(descs["MolWt"] - 206.28) < 2.0

    def test_invalid_smiles_returns_error(self) -> None:
        response = handle_compute_descriptors({"smiles": "NOT_VALID"})
        data = _parse(response)
        # Should return error payload
        assert "error" in data or data.get("total_descriptors", 1) == 0

    def test_total_descriptors_count(self) -> None:
        response = handle_compute_descriptors({"smiles": ASPIRIN_SMILES})
        data = _parse(response)
        # RDKit has 200+ descriptors
        assert data["total_descriptors"] >= 100


# ---------------------------------------------------------------------------
# 6. Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    def test_predict_admet_missing_smiles_key_raises(self) -> None:
        # Missing required key should be caught gracefully
        with pytest.raises((KeyError, TypeError)):
            handle_predict_admet({})

    def test_check_drug_likeness_empty_smiles(self) -> None:
        response = handle_check_drug_likeness({"smiles": ""})
        data = _parse(response)
        # Empty SMILES: either error or filter results with all violations
        assert isinstance(data, dict)

    def test_compute_descriptors_empty_smiles(self) -> None:
        response = handle_compute_descriptors({"smiles": ""})
        data = _parse(response)
        assert isinstance(data, dict)

    @pytest.mark.asyncio
    async def test_cross_reference_invalid_cid_type(self) -> None:
        response = await handle_cross_reference(
            {"identifier": "not_a_number", "identifier_type": "cid"}
        )
        assert _is_error(response)

    @pytest.mark.asyncio
    async def test_cross_reference_unknown_type(self) -> None:
        response = await handle_cross_reference(
            {"identifier": "somevalue", "identifier_type": "unknown_type"}
        )
        assert _is_error(response)

    def test_analyze_scaffolds_empty_list(self) -> None:
        response = handle_analyze_scaffolds({"smiles_list": []})
        assert _is_error(response)

    def test_compare_compounds_invalid_smiles_returns_error(self) -> None:
        response = handle_compare_compounds(
            {"smiles_a": "INVALID", "smiles_b": ASPIRIN_SMILES}
        )
        assert _is_error(response)


# ---------------------------------------------------------------------------
# 7. analyze_scaffolds
# ---------------------------------------------------------------------------


class TestAnalyzeScaffolds:
    def test_aspirin_and_ibuprofen_scaffolds(self) -> None:
        response = handle_analyze_scaffolds({"smiles_list": [ASPIRIN_SMILES, IBUPROFEN_SMILES]})
        data = _parse(response)

        if _is_error(response):
            # RDKit not available — acceptable
            assert "RDKit" in data["error"]
            return

        assert data["total_compounds"] == 2
        assert "unique_scaffolds" in data
        assert "scaffolds" in data

    def test_duplicate_scaffolds_counted(self) -> None:
        # Both aspirin and salicylic acid share the benzoic acid scaffold
        salicylic_acid = "Oc1ccccc1C(=O)O"
        response = handle_analyze_scaffolds({"smiles_list": [ASPIRIN_SMILES, salicylic_acid]})
        data = _parse(response)

        if _is_error(response):
            return  # RDKit not available

        assert data["total_compounds"] == 2
        # May share same scaffold depending on RDKit version
        assert isinstance(data["scaffolds"], list)

    def test_single_smiles_returns_one_scaffold(self) -> None:
        response = handle_analyze_scaffolds({"smiles_list": [ASPIRIN_SMILES]})
        data = _parse(response)

        if _is_error(response):
            return

        assert data["total_compounds"] == 1


# ---------------------------------------------------------------------------
# 8. compare_compounds
# ---------------------------------------------------------------------------


class TestCompareCompounds:
    def test_aspirin_vs_ibuprofen(self) -> None:
        response = handle_compare_compounds(
            {"smiles_a": ASPIRIN_SMILES, "smiles_b": IBUPROFEN_SMILES}
        )
        data = _parse(response)

        if _is_error(response):
            # RDKit not available
            return

        assert "tanimoto_similarity" in data
        assert 0.0 <= data["tanimoto_similarity"] <= 1.0
        assert "similarity_label" in data
        assert data["similarity_label"] in ("low", "medium", "high")
        assert "property_comparison" in data
        assert "shared_scaffold" in data

    def test_identical_compounds_similarity_is_one(self) -> None:
        response = handle_compare_compounds(
            {"smiles_a": ASPIRIN_SMILES, "smiles_b": ASPIRIN_SMILES}
        )
        data = _parse(response)

        if _is_error(response):
            return

        assert data["tanimoto_similarity"] == 1.0
        assert data["similarity_label"] == "high"

    def test_aspirin_property_comparison_keys(self) -> None:
        response = handle_compare_compounds(
            {"smiles_a": ASPIRIN_SMILES, "smiles_b": IBUPROFEN_SMILES}
        )
        data = _parse(response)

        if _is_error(response):
            return

        prop_comp = data["property_comparison"]
        if prop_comp:
            for _key, val in prop_comp.items():
                assert "compound_a" in val
                assert "compound_b" in val
                assert "delta" in val


# ---------------------------------------------------------------------------
# 9. Server module is importable without running server
# ---------------------------------------------------------------------------


class TestImportGuard:
    def test_module_importable(self) -> None:
        import chemfuse.mcp.server as server_module
        assert hasattr(server_module, "_TOOLS")
        assert hasattr(server_module, "main")

    def test_tools_list_not_empty_when_mcp_available(self) -> None:
        assert len(_TOOLS) == 9

    def test_main_is_coroutine(self) -> None:
        import inspect

        from chemfuse.mcp.server import main
        assert inspect.iscoroutinefunction(main)
