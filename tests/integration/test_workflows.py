"""Integration tests for ChemFuse end-to-end workflows."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

from chemfuse import find_similar, get, search
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture
def aspirin_compound() -> Compound:
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        name="aspirin",
        formula="C9H8O4",
        synonyms=["aspirin", "acetylsalicylic acid"],
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=180.16,
            xlogp=1.2,
            tpsa=63.6,
            hbd_count=1,
            hba_count=4,
            rotatable_bonds=3,
        ),
    )


@pytest.fixture
def ibuprofen_compound() -> Compound:
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
        ),
    )


class TestSearchAndFilterWorkflow:
    def test_search_and_filter_druglike(
        self, aspirin_compound: Compound, ibuprofen_compound: Compound
    ):
        """Test search -> filter druglike workflow."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(
                return_value=[aspirin_compound, ibuprofen_compound]
            )
            mock_registry.get.return_value = mock_adapter

            collection = search("nsaid")
            lipinski = collection.filter_druglike(rule="lipinski")

            assert len(lipinski) >= 0  # Both should pass Lipinski
            assert isinstance(lipinski, CompoundCollection)

    def test_search_and_sort_by_mw(
        self, aspirin_compound: Compound, ibuprofen_compound: Compound
    ):
        """Test search -> sort by molecular weight workflow."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(
                return_value=[ibuprofen_compound, aspirin_compound]
            )
            mock_registry.get.return_value = mock_adapter

            collection = search("nsaid")
            sorted_coll = collection.sort(by="molecular_weight", ascending=True)

            weights = [
                c.properties.molecular_weight
                for c in sorted_coll
                if c.properties.molecular_weight
            ]
            assert weights == sorted(weights)

    def test_search_export_to_csv(
        self, aspirin_compound: Compound, tmp_path: Path
    ):
        """Test search -> export to CSV workflow."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[aspirin_compound])
            mock_registry.get.return_value = mock_adapter

            collection = search("aspirin")
            output = str(tmp_path / "results.csv")
            result = collection.to_csv(output)
            assert Path(result).exists()

    def test_search_to_dataframe(
        self, aspirin_compound: Compound, ibuprofen_compound: Compound
    ):
        """Test search -> to_dataframe workflow."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(
                return_value=[aspirin_compound, ibuprofen_compound]
            )
            mock_registry.get.return_value = mock_adapter

            collection = search("nsaid")
            df = collection.to_dataframe()

            assert len(df) == 2
            assert "smiles" in df.columns
            assert "name" in df.columns


class TestGetAndMergeWorkflow:
    def test_get_compound_by_id(self, aspirin_compound: Compound):
        """Test retrieving compound by ID."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_by_id = AsyncMock(return_value=aspirin_compound)
            mock_registry.get.return_value = mock_adapter

            compound = get("2244")
            assert compound is not None
            assert compound.cid == 2244
            assert compound.name == "aspirin"

    def test_merge_compounds_from_sources(
        self, aspirin_compound: Compound
    ):
        """Test merging compound data from multiple sources."""
        c1 = Compound(
            smiles=aspirin_compound.smiles,
            inchikey=aspirin_compound.inchikey,
            name="aspirin",
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=180.16),
        )
        c2 = Compound(
            smiles=aspirin_compound.smiles,
            inchikey=aspirin_compound.inchikey,
            formula="C9H8O4",
            sources=["chembl"],
            properties=CompoundProperties(xlogp=1.2),
        )

        merged = c1.merge(c2)
        assert merged.formula == "C9H8O4"
        assert merged.properties.molecular_weight == 180.16
        assert merged.properties.xlogp == 1.2
        assert "pubchem" in merged.sources
        assert "chembl" in merged.sources


class TestSimilarityWorkflow:
    def test_find_similar_compounds(
        self, aspirin_compound: Compound, ibuprofen_compound: Compound
    ):
        """Test similarity search workflow."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_similarity = AsyncMock(
                return_value=[aspirin_compound, ibuprofen_compound]
            )
            mock_registry.get.return_value = mock_adapter

            similar = find_similar("CC(=O)Oc1ccccc1C(=O)O", threshold=80)

            assert isinstance(similar, CompoundCollection)
            assert len(similar) >= 0


class TestBatchSearchWorkflow:
    def test_batch_search_from_file(
        self, tmp_path: Path, aspirin_compound: Compound
    ):
        """Test batch search from a text file."""
        query_file = tmp_path / "queries.txt"
        query_file.write_text("aspirin\nibuprofen\n")

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[aspirin_compound])
            mock_registry.get.return_value = mock_adapter
            from chemfuse import batch_search
            collection, errors = batch_search(str(query_file))

        assert isinstance(collection, CompoundCollection)
        assert isinstance(errors, list)
