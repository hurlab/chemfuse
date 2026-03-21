"""Tests for chemfuse.sources.pubchem (PubChemAdapter)."""

from __future__ import annotations

import httpx
import pytest
import respx

from chemfuse.core.exceptions import NotFoundError, TimeoutError
from chemfuse.sources.pubchem import PUBCHEM_BASE_URL, PubChemAdapter


@pytest.fixture
def adapter() -> PubChemAdapter:
    """PubChemAdapter with high rate limit for fast tests."""
    return PubChemAdapter(timeout=5.0)


@pytest.fixture
def property_table_response() -> dict:
    return {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 2244,
                    "MolecularFormula": "C9H8O4",
                    "MolecularWeight": 180.16,
                    "ExactMass": 180.042,
                    "CanonicalSMILES": "CC(=O)Oc1ccccc1C(=O)O",
                    "IsomericSMILES": "CC(=O)Oc1ccccc1C(=O)O",
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


@pytest.fixture
def synonyms_response() -> dict:
    return {
        "InformationList": {
            "Information": [
                {
                    "CID": 2244,
                    "Synonym": ["aspirin", "acetylsalicylic acid", "50-78-2"],
                }
            ]
        }
    }


@pytest.fixture
def identifier_list_response() -> dict:
    return {"IdentifierList": {"CID": [2244, 3672, 2519]}}


class TestPubChemAdapterInit:
    def test_instantiation(self, adapter: PubChemAdapter):
        assert adapter is not None
        assert adapter.name == "pubchem"
        assert adapter.base_url == PUBCHEM_BASE_URL

    def test_is_available(self, adapter: PubChemAdapter):
        assert adapter.is_available() is True


class TestSearchByName:
    @respx.mock
    async def test_search_by_name_returns_compounds(
        self, adapter: PubChemAdapter, property_table_response: dict, synonyms_response: dict
    ):
        respx.get(url__regex=r".*/name/aspirin/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )
        respx.get(url__regex=r".*/synonyms/.*").mock(
            return_value=httpx.Response(200, json=synonyms_response)
        )

        compounds = await adapter.search("aspirin", query_type="name")
        assert len(compounds) == 1
        assert compounds[0].cid == 2244
        assert compounds[0].name is not None

    @respx.mock
    async def test_search_by_name_not_found(self, adapter: PubChemAdapter):
        respx.get(url__regex=r".*/name/.*").mock(
            return_value=httpx.Response(
                404, json={"Fault": {"Message": "No CID found"}}
            )
        )
        with pytest.raises(NotFoundError):
            await adapter.search("nonexistent_compound_xyz", query_type="name")

    @respx.mock
    async def test_search_uses_default_type(
        self, adapter: PubChemAdapter, property_table_response: dict, synonyms_response: dict
    ):
        respx.get(url__regex=r".*/name/aspirin/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )
        respx.get(url__regex=r".*/synonyms/.*").mock(
            return_value=httpx.Response(200, json=synonyms_response)
        )
        # Default query_type should be "name"
        compounds = await adapter.search("aspirin")
        assert len(compounds) >= 1


class TestSearchByCID:
    @respx.mock
    async def test_search_by_cid_returns_compound(
        self, adapter: PubChemAdapter, property_table_response: dict, synonyms_response: dict
    ):
        respx.get(url__regex=r".*/cid/2244/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )
        respx.get(url__regex=r".*/cid/2244/synonyms/.*").mock(
            return_value=httpx.Response(200, json=synonyms_response)
        )

        compounds = await adapter.search("2244", query_type="cid")
        assert len(compounds) == 1
        assert compounds[0].cid == 2244

    @respx.mock
    async def test_search_by_identifier_with_digit_string(
        self, adapter: PubChemAdapter, property_table_response: dict, synonyms_response: dict
    ):
        respx.get(url__regex=r".*/cid/2244/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )
        respx.get(url__regex=r".*/synonyms/.*").mock(
            return_value=httpx.Response(200, json=synonyms_response)
        )
        compounds = await adapter.search("2244", query_type="identifier")
        assert len(compounds) == 1


class TestSearchBySMILES:
    @respx.mock
    async def test_search_by_smiles(
        self, adapter: PubChemAdapter, property_table_response: dict
    ):
        respx.get(url__regex=r".*/smiles/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.search("CC(=O)Oc1ccccc1C(=O)O", query_type="smiles")
        assert len(compounds) == 1
        assert compounds[0].smiles == "CC(=O)Oc1ccccc1C(=O)O"


class TestSearchByFormula:
    @respx.mock
    async def test_search_by_formula(
        self, adapter: PubChemAdapter, property_table_response: dict, identifier_list_response: dict
    ):
        respx.get(url__regex=r".*/fastformula/.*").mock(
            return_value=httpx.Response(200, json=identifier_list_response)
        )
        respx.get(url__regex=r".*/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.search("C9H8O4", query_type="formula")
        assert isinstance(compounds, list)


class TestSearchByInChI:
    @respx.mock
    async def test_search_by_inchi(
        self, adapter: PubChemAdapter, property_table_response: dict
    ):
        respx.post(url__regex=r".*/inchi/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.search(
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            query_type="inchi",
        )
        assert isinstance(compounds, list)


class TestGetByID:
    @respx.mock
    async def test_get_by_id_returns_compound(
        self, adapter: PubChemAdapter, property_table_response: dict, synonyms_response: dict
    ):
        respx.get(url__regex=r".*/cid/2244/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )
        respx.get(url__regex=r".*/synonyms/.*").mock(
            return_value=httpx.Response(200, json=synonyms_response)
        )

        compound = await adapter.get_by_id("2244")
        assert compound is not None
        assert compound.cid == 2244

    @respx.mock
    async def test_get_by_id_not_found_returns_none(self, adapter: PubChemAdapter):
        respx.get(url__regex=r".*/cid/99999/.*").mock(
            return_value=httpx.Response(404, json={"Fault": {"Message": "Not found"}})
        )
        compound = await adapter.get_by_id("99999")
        assert compound is None


class TestGetProperties:
    @respx.mock
    async def test_get_properties_returns_dict(
        self, adapter: PubChemAdapter, property_table_response: dict
    ):
        respx.get(url__regex=r".*/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        props = await adapter.get_properties("2244")
        assert isinstance(props, dict)
        assert props.get("CID") == 2244 or props.get("MolecularWeight") is not None


class TestGetSimilarity:
    @respx.mock
    async def test_similarity_search_direct_list(
        self, adapter: PubChemAdapter, property_table_response: dict, identifier_list_response: dict
    ):
        respx.get(url__regex=r".*/similarity/smiles/.*").mock(
            return_value=httpx.Response(200, json=identifier_list_response)
        )
        respx.get(url__regex=r".*/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.get_similarity("CC(=O)Oc1ccccc1C(=O)O", threshold=90)
        assert isinstance(compounds, list)

    @respx.mock
    async def test_similarity_search_with_list_key(
        self,
        adapter: PubChemAdapter,
        property_table_response: dict,
        identifier_list_response: dict,
    ):
        waiting_response = {"Waiting": {"ListKey": "test-key-123"}}

        respx.get(url__regex=r".*/similarity/smiles/.*").mock(
            return_value=httpx.Response(200, json=waiting_response)
        )
        respx.get(url__regex=r".*/listkey/.*").mock(
            return_value=httpx.Response(200, json=identifier_list_response)
        )
        respx.get(url__regex=r".*/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.get_similarity("CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(compounds, list)

    @respx.mock
    async def test_similarity_source_error_returns_empty(self, adapter: PubChemAdapter):
        respx.get(url__regex=r".*/similarity/smiles/.*").mock(
            return_value=httpx.Response(500, json={"error": "server error"})
        )
        compounds = await adapter.get_similarity("INVALID_SMILES")
        assert compounds == []


class TestGetSubstructure:
    @respx.mock
    async def test_substructure_search(
        self, adapter: PubChemAdapter, property_table_response: dict, identifier_list_response: dict
    ):
        respx.get(url__regex=r".*/substructure/smiles/.*").mock(
            return_value=httpx.Response(200, json=identifier_list_response)
        )
        respx.get(url__regex=r".*/property/.*").mock(
            return_value=httpx.Response(200, json=property_table_response)
        )

        compounds = await adapter.get_substructure("c1ccccc1")
        assert isinstance(compounds, list)


class TestGetBioassays:
    @respx.mock
    async def test_get_bioassays_returns_list(self, adapter: PubChemAdapter):
        assay_response = {
            "Table": {
                "Columns": {"Column": ["AID", "Name", "Outcome"]},
                "Row": [
                    {"Cell": ["12345", "Test Assay", "Active"]},
                ],
            }
        }
        respx.get(url__regex=r".*/assaysummary/.*").mock(
            return_value=httpx.Response(200, json=assay_response)
        )

        assays = await adapter.get_bioassays(2244)
        assert isinstance(assays, list)
        assert len(assays) == 1
        assert assays[0]["AID"] == "12345"

    @respx.mock
    async def test_get_bioassays_not_found_returns_empty(self, adapter: PubChemAdapter):
        respx.get(url__regex=r".*/assaysummary/.*").mock(
            return_value=httpx.Response(404, json={"Fault": {"Message": "Not found"}})
        )
        assays = await adapter.get_bioassays(99999)
        assert assays == []


class TestParsePropDict:
    def test_parse_complete_prop_dict(self, adapter: PubChemAdapter):
        prop = {
            "CID": 2244,
            "MolecularFormula": "C9H8O4",
            "MolecularWeight": 180.16,
            "ExactMass": 180.042,
            "CanonicalSMILES": "CC(=O)Oc1ccccc1C(=O)O",
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
        compound = adapter._prop_dict_to_compound(prop)
        assert compound.cid == 2244
        assert compound.formula == "C9H8O4"
        assert compound.smiles == "CC(=O)Oc1ccccc1C(=O)O"
        assert compound.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        assert compound.properties.molecular_weight == 180.16
        assert compound.properties.xlogp == 1.2
        assert compound.properties.hbd_count == 1
        assert "pubchem" in compound.sources

    def test_parse_minimal_prop_dict(self, adapter: PubChemAdapter):
        prop = {"CID": 100, "CanonicalSMILES": "C"}
        compound = adapter._prop_dict_to_compound(prop)
        assert compound.cid == 100
        assert compound.smiles == "C"

    def test_falls_back_to_isomeric_smiles(self, adapter: PubChemAdapter):
        prop = {"CID": 100, "IsomericSMILES": "[C@@H](F)(Cl)Br"}
        compound = adapter._prop_dict_to_compound(prop)
        assert compound.smiles == "[C@@H](F)(Cl)Br"


class TestResolveListKeyOrIds:
    async def test_resolves_direct_identifier_list(self, adapter: PubChemAdapter):
        data = {"IdentifierList": {"CID": [1, 2, 3]}}
        cids = await adapter._resolve_list_key_or_ids(data)
        assert cids == [1, 2, 3]

    async def test_returns_empty_for_unknown_format(self, adapter: PubChemAdapter):
        data = {"Unknown": "data"}
        cids = await adapter._resolve_list_key_or_ids(data)
        assert cids == []

    async def test_returns_empty_for_non_dict(self, adapter: PubChemAdapter):
        cids = await adapter._resolve_list_key_or_ids([1, 2, 3])
        assert cids == []


class TestPollListKey:
    @respx.mock
    async def test_poll_list_key_success(
        self, adapter: PubChemAdapter, identifier_list_response: dict
    ):
        respx.get(url__regex=r".*/listkey/.*").mock(
            return_value=httpx.Response(200, json=identifier_list_response)
        )
        cids = await adapter._poll_list_key("test-key-123", max_wait=10.0, poll_interval=0.01)
        assert cids == [2244, 3672, 2519]

    @respx.mock
    async def test_poll_list_key_timeout(self, adapter: PubChemAdapter):
        waiting = {"Waiting": {"ListKey": "test-key-123"}}
        respx.get(url__regex=r".*/listkey/.*").mock(
            return_value=httpx.Response(200, json=waiting)
        )
        with pytest.raises(TimeoutError):
            await adapter._poll_list_key("test-key-123", max_wait=0.1, poll_interval=0.01)
