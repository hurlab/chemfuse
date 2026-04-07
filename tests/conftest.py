"""Shared pytest fixtures for ChemFuse tests."""

from __future__ import annotations

import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture(autouse=True)
def _reset_registry_adapters():
    """Reset cached adapter instances between tests.

    The registry injects a shared Cache into adapters. Cached API responses
    from previous tests (or real calls) can bypass respx mocks. Clearing
    cached instances and using a disabled cache ensures test isolation.
    """
    from chemfuse.core.cache import Cache
    from chemfuse.sources import registry
    registry._adapters.clear()
    registry._cache_instance = Cache(enabled=False)
    yield
    registry._adapters.clear()
    registry._cache_instance = None


@pytest.fixture
def aspirin_props() -> CompoundProperties:
    """CompoundProperties for aspirin."""
    return CompoundProperties(
        molecular_weight=180.16,
        exact_mass=180.0423,
        xlogp=1.2,
        tpsa=63.6,
        hbd_count=1,
        hba_count=4,
        rotatable_bonds=3,
        heavy_atom_count=13,
        complexity=212.0,
    )


@pytest.fixture
def aspirin(aspirin_props: CompoundProperties) -> Compound:
    """Compound object for aspirin."""
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchi="InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        name="aspirin",
        formula="C9H8O4",
        synonyms=["aspirin", "acetylsalicylic acid", "2-acetoxybenzoic acid"],
        sources=["pubchem"],
        properties=aspirin_props,
    )


@pytest.fixture
def ibuprofen() -> Compound:
    """Compound object for ibuprofen."""
    return Compound(
        cid=3672,
        smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        inchikey="HEFNNWSXXWATRW-UHFFFAOYSA-N",
        name="ibuprofen",
        formula="C13H18O2",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=206.28,
            xlogp=3.5,
            tpsa=37.3,
            hbd_count=1,
            hba_count=2,
            rotatable_bonds=4,
            heavy_atom_count=15,
        ),
    )


@pytest.fixture
def caffeine() -> Compound:
    """Compound object for caffeine."""
    return Compound(
        cid=2519,
        smiles="Cn1c(=O)c2c(ncn2C)n(c1=O)C",
        inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        name="caffeine",
        formula="C8H10N4O2",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=194.19,
            xlogp=-0.07,
            tpsa=58.4,
            hbd_count=0,
            hba_count=6,
            rotatable_bonds=0,
            heavy_atom_count=14,
        ),
    )


@pytest.fixture
def sample_collection(aspirin: Compound, ibuprofen: Compound, caffeine: Compound) -> CompoundCollection:
    """A collection with aspirin, ibuprofen, and caffeine."""
    return CompoundCollection(
        compounds=[aspirin, ibuprofen, caffeine],
        query="test",
        sources=["pubchem"],
    )


@pytest.fixture
def empty_collection() -> CompoundCollection:
    """Empty CompoundCollection."""
    return CompoundCollection()


# PubChem API response fixtures

@pytest.fixture
def pubchem_property_table_aspirin() -> dict:
    """Mock PubChem PropertyTable response for aspirin."""
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
                    "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
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
def pubchem_synonyms_aspirin() -> dict:
    """Mock PubChem synonyms response for aspirin."""
    return {
        "InformationList": {
            "Information": [
                {
                    "CID": 2244,
                    "Synonym": [
                        "aspirin",
                        "acetylsalicylic acid",
                        "2-acetoxybenzoic acid",
                        "50-78-2",
                    ],
                }
            ]
        }
    }


@pytest.fixture
def pubchem_identifier_list() -> dict:
    """Mock PubChem IdentifierList response."""
    return {
        "IdentifierList": {
            "CID": [2244, 3672, 2519],
        }
    }


@pytest.fixture
def pubchem_waiting_response() -> dict:
    """Mock PubChem ListKey waiting response."""
    return {
        "Waiting": {
            "ListKey": "test-list-key-12345",
        }
    }
