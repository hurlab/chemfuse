"""Tests for BindingDB source adapter."""

from __future__ import annotations

import pytest
import respx
from httpx import Response

from chemfuse.models.bioactivity import BindingMeasurement
from chemfuse.sources.bindingdb import BindingDBAdapter, _parse_affinity

BINDINGDB_BASE = "https://www.bindingdb.org/axis2/services/BDBService"

# Mock JSON responses
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(O)=O"

BINDINGDB_JSON_RESPONSE = {
    "affinities": [
        {
            "target": "Cyclooxygenase-1",
            "uniprot": "P23219",
            "ki": "1670",
            "kd": None,
            "ic50": None,
            "ec50": None,
        },
        {
            "target": "Cyclooxygenase-2",
            "uniprot": "P35354",
            "ki": None,
            "kd": None,
            "ic50": "278000",
            "ec50": None,
        },
    ]
}

BINDINGDB_JSON_INEQUALITY = {
    "affinities": [
        {
            "target": "Target A",
            "ki": ">10000",
            "ic50": "<100",
        }
    ]
}

BINDINGDB_XML_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<response>
  <affinity>
    <target>Cyclooxygenase-1</target>
    <uniprot>P23219</uniprot>
    <ki>1670</ki>
    <ic50></ic50>
  </affinity>
  <affinity>
    <target>Cyclooxygenase-2</target>
    <uniprot>P35354</uniprot>
    <ki></ki>
    <ic50>278000</ic50>
  </affinity>
</response>"""

COX1_UNIPROT = "P23219"
UNIPROT_RESPONSE = {
    "affinities": [
        {
            "target": "Cyclooxygenase-1",
            "uniprot": "P23219",
            "ki": "500",
        }
    ]
}


class TestParseAffinity:
    """Tests for the _parse_affinity helper function."""

    def test_plain_numeric(self) -> None:
        """Plain numeric string is parsed to float with no relation."""
        value, relation = _parse_affinity("1670")
        assert value == pytest.approx(1670.0)
        assert relation is None

    def test_float_string(self) -> None:
        """Decimal numeric string is parsed correctly."""
        value, relation = _parse_affinity("1.67")
        assert value == pytest.approx(1.67)
        assert relation is None

    def test_greater_than(self) -> None:
        """'>10000' is parsed to value=10000 and relation='>'."""
        value, relation = _parse_affinity(">10000")
        assert value == pytest.approx(10000.0)
        assert relation == ">"

    def test_less_than(self) -> None:
        """'<100' is parsed to value=100 and relation='<'."""
        value, relation = _parse_affinity("<100")
        assert value == pytest.approx(100.0)
        assert relation == "<"

    def test_greater_than_equal(self) -> None:
        """'>=500' is parsed correctly."""
        value, relation = _parse_affinity(">=500")
        assert value == pytest.approx(500.0)
        assert relation == ">="

    def test_tilde_prefix(self) -> None:
        """'~100' is parsed with relation='~'."""
        value, relation = _parse_affinity("~100")
        assert value == pytest.approx(100.0)
        assert relation == "~"

    def test_none_returns_none(self) -> None:
        """None input returns (None, None)."""
        value, relation = _parse_affinity(None)
        assert value is None
        assert relation is None

    def test_empty_string_returns_none(self) -> None:
        """Empty string returns (None, None)."""
        value, relation = _parse_affinity("")
        assert value is None
        assert relation is None

    def test_na_string_returns_none(self) -> None:
        """'N/A' strings return (None, None)."""
        for na in ("na", "N/A", "none", "-"):
            value, relation = _parse_affinity(na)
            assert value is None
            assert relation is None

    def test_numeric_int(self) -> None:
        """Numeric int input is parsed to float."""
        value, relation = _parse_affinity(500)
        assert value == pytest.approx(500.0)
        assert relation is None


class TestBindingDBSearchBySmiles:
    """Tests for BindingDB search_by_smiles."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles_returns_measurements(self) -> None:
        """search_by_smiles returns BindingMeasurement objects."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(200, json=BINDINGDB_JSON_RESPONSE)
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert len(measurements) == 2
        assert all(isinstance(m, BindingMeasurement) for m in measurements)

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles_correct_values(self) -> None:
        """search_by_smiles parses Ki and IC50 values correctly."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(200, json=BINDINGDB_JSON_RESPONSE)
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles(ASPIRIN_SMILES)
        # First measurement: COX-1 with Ki=1670
        cox1 = next(m for m in measurements if m.target_name == "Cyclooxygenase-1")
        assert cox1.ki == pytest.approx(1670.0)
        assert cox1.target_uniprot == "P23219"
        # Second measurement: COX-2 with IC50=278000
        cox2 = next(m for m in measurements if m.target_name == "Cyclooxygenase-2")
        assert cox2.ic50 == pytest.approx(278000.0)

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles_inequality_values(self) -> None:
        """search_by_smiles parses inequality affinity values."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(200, json=BINDINGDB_JSON_INEQUALITY)
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert len(measurements) == 1
        m = measurements[0]
        assert m.ki == pytest.approx(10000.0)
        assert m.ki_relation == ">"
        assert m.ic50 == pytest.approx(100.0)
        assert m.ic50_relation == "<"

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles_empty_result(self) -> None:
        """search_by_smiles returns empty list when no data found."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(200, json={"affinities": []})
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles("CC")
        assert measurements == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_smiles_api_error_returns_empty(self) -> None:
        """search_by_smiles returns empty list on API error."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(500, text="Server Error")
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert measurements == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_source_is_set_correctly(self) -> None:
        """Measurements have source='bindingdb'."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(200, json=BINDINGDB_JSON_RESPONSE)
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_smiles(ASPIRIN_SMILES)
        assert all(m.source == "bindingdb" for m in measurements)


class TestBindingDBSearchByTarget:
    """Tests for BindingDB search_by_target."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_target_returns_measurements(self) -> None:
        """search_by_target returns BindingMeasurement objects."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsByUniProt").mock(
            return_value=Response(200, json=UNIPROT_RESPONSE)
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_target(COX1_UNIPROT)
        assert len(measurements) == 1
        assert measurements[0].ki == pytest.approx(500.0)

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_target_empty_result(self) -> None:
        """search_by_target returns empty list when no data found."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsByUniProt").mock(
            return_value=Response(200, json={"affinities": []})
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_target("Q99999")
        assert measurements == []

    @pytest.mark.asyncio
    @respx.mock
    async def test_search_by_target_error_returns_empty(self) -> None:
        """search_by_target returns empty list on API failure."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsByUniProt").mock(
            return_value=Response(503, text="Service Unavailable")
        )
        adapter = BindingDBAdapter()
        measurements = await adapter.search_by_target("P23219")
        assert measurements == []


class TestBindingDBXMLParsing:
    """Tests for BindingDB XML response parsing."""

    @pytest.mark.asyncio
    @respx.mock
    async def test_xml_response_is_parsed(self) -> None:
        """XML response is parsed into BindingMeasurement objects."""
        respx.get(f"{BINDINGDB_BASE}/getLigandsBySmiles").mock(
            return_value=Response(
                200,
                content=BINDINGDB_XML_RESPONSE.encode(),
                headers={"content-type": "text/xml"},
            )
        )
        adapter = BindingDBAdapter()
        # The _parse_response handles XML when data is a string
        measurements = adapter._parse_xml(BINDINGDB_XML_RESPONSE)
        assert len(measurements) == 2
        cox1 = next(m for m in measurements if m.target_name == "Cyclooxygenase-1")
        assert cox1.ki == pytest.approx(1670.0)

    def test_malformed_xml_returns_empty(self) -> None:
        """Malformed XML returns empty list."""
        adapter = BindingDBAdapter()
        result = adapter._parse_xml("<not valid xml")
        assert result == []


class TestBindingDBAdapterMisc:
    """Miscellaneous BindingDB adapter tests."""

    def test_adapter_name(self) -> None:
        """Adapter has correct name."""
        adapter = BindingDBAdapter()
        assert adapter.name == "bindingdb"

    def test_is_available(self) -> None:
        """is_available always returns True."""
        adapter = BindingDBAdapter()
        assert adapter.is_available() is True

    @pytest.mark.asyncio
    async def test_search_returns_empty(self) -> None:
        """search always returns empty list."""
        adapter = BindingDBAdapter()
        result = await adapter.search("CC")
        assert result == []

    @pytest.mark.asyncio
    async def test_get_by_id_returns_none(self) -> None:
        """get_by_id always returns None."""
        adapter = BindingDBAdapter()
        result = await adapter.get_by_id("12345")
        assert result is None

    @pytest.mark.asyncio
    async def test_get_properties_returns_empty(self) -> None:
        """get_properties always returns empty dict."""
        adapter = BindingDBAdapter()
        result = await adapter.get_properties("P23219")
        assert result == {}

    @pytest.mark.asyncio
    @respx.mock
    async def test_rate_limiting_attr(self) -> None:
        """Adapter has rate_limit attribute."""
        adapter = BindingDBAdapter()
        assert adapter.rate_limit == 1.0


class TestBindingDBParseEdgeCases:
    """Edge case tests for BindingDB response parsing."""

    def test_parse_response_with_list_data(self) -> None:
        """_parse_response handles list input directly."""
        adapter = BindingDBAdapter()
        items = [
            {"target": "COX-1", "uniprot": "P23219", "ki": "500"},
        ]
        result = adapter._parse_response(items)
        assert len(result) == 1
        assert result[0].target_name == "COX-1"

    def test_parse_response_with_string_triggers_xml(self) -> None:
        """_parse_response with string input falls through to XML parsing."""
        adapter = BindingDBAdapter()
        xml = """<response><affinity><target>COX-1</target><ki>500</ki></affinity></response>"""
        result = adapter._parse_response(xml)
        assert len(result) == 1
        assert result[0].ki == pytest.approx(500.0)

    def test_parse_json_nested_affinities_key(self) -> None:
        """_parse_json handles nested dict with affinities key."""
        adapter = BindingDBAdapter()
        data = {
            "getLigandsBySmilesResponse": {
                "affinities": [
                    {"target": "COX-1", "ki": "1000"},
                ]
            }
        }
        result = adapter._parse_json(data)
        assert len(result) == 1
        assert result[0].ki == pytest.approx(1000.0)

    def test_parse_json_non_list_measurements_raw(self) -> None:
        """_parse_json returns empty list when affinities value is not a list."""
        adapter = BindingDBAdapter()
        data = {"affinities": "not_a_list"}
        result = adapter._parse_json(data)
        assert result == []

    def test_parse_json_list_skips_non_dict_items(self) -> None:
        """_parse_json_list skips items that are not dicts."""
        adapter = BindingDBAdapter()
        items = [
            "not_a_dict",
            {"target": "COX-1", "ki": "500"},
        ]
        result = adapter._parse_json_list(items)
        assert len(result) == 1

    def test_parse_measurement_dict_skips_empty_record(self) -> None:
        """_parse_measurement_dict returns None for item with no target and no values."""
        adapter = BindingDBAdapter()
        result = adapter._parse_measurement_dict({})
        assert result is None

    def test_xml_empty_affinity_element_skipped(self) -> None:
        """XML affinity with no target and no values is skipped."""
        adapter = BindingDBAdapter()
        xml = """<response>
            <affinity></affinity>
            <affinity><target>COX-1</target><ki>500</ki></affinity>
        </response>"""
        result = adapter._parse_xml(xml)
        assert len(result) == 1
        assert result[0].target_name == "COX-1"

    def test_parse_affinity_fallback_float_parse(self) -> None:
        """_parse_affinity handles strings that don't match the inequality regex."""
        # A string like "1e3" (scientific notation) won't match the regex
        # and falls back to float() which succeeds
        value, relation = _parse_affinity("1e3")
        assert value == pytest.approx(1000.0)
        assert relation is None

    def test_parse_affinity_non_numeric_fallback_returns_none(self) -> None:
        """_parse_affinity returns None for completely non-numeric strings."""
        value, relation = _parse_affinity("not_a_number")
        assert value is None
        assert relation is None
