"""Tests for the ChemFuse public API (chemfuse.__init__)."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest

import chemfuse
from chemfuse import (
    Compound,
    CompoundCollection,
    CompoundProperties,
    __version__,
    find_similar,
    find_similar_async,
    get,
    get_async,
    search,
    search_async,
)
from chemfuse.models.collection import CompoundCollection as CollectionType


@pytest.fixture
def aspirin_compound() -> Compound:
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem"],
        properties=CompoundProperties(molecular_weight=180.16),
    )


@pytest.fixture
def aspirin_collection(aspirin_compound: Compound) -> CompoundCollection:
    return CompoundCollection(
        compounds=[aspirin_compound],
        query="aspirin",
        sources=["pubchem"],
    )


class TestVersion:
    def test_version_is_string(self):
        assert isinstance(__version__, str)

    def test_version_format(self):
        parts = __version__.split(".")
        assert len(parts) >= 2


class TestPublicAPIExports:
    def test_compound_importable(self):
        assert chemfuse.Compound is not None

    def test_compound_collection_importable(self):
        assert chemfuse.CompoundCollection is not None

    def test_compound_properties_importable(self):
        assert chemfuse.CompoundProperties is not None

    def test_search_importable(self):
        assert callable(chemfuse.search)

    def test_get_importable(self):
        assert callable(chemfuse.get)

    def test_find_similar_importable(self):
        assert callable(chemfuse.find_similar)

    def test_cross_reference_importable(self):
        assert callable(chemfuse.cross_reference)

    def test_batch_search_importable(self):
        assert callable(chemfuse.batch_search)


class TestSearchAsync:
    async def test_search_async_returns_collection(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=aspirin_collection.compounds)
            mock_registry.get.return_value = mock_adapter

            result = await search_async("aspirin")
            assert isinstance(result, CollectionType)

    async def test_search_async_default_source_is_pubchem(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            await search_async("aspirin")
            mock_registry.get.assert_called_with("pubchem")

    async def test_search_async_multiple_sources(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=aspirin_collection.compounds)
            mock_registry.get.return_value = mock_adapter

            result = await search_async("aspirin", sources=["pubchem"])
            assert isinstance(result, CollectionType)

    async def test_search_async_source_error_returns_partial(self):
        """A source error in one source should not fail the whole search."""
        from chemfuse.core.exceptions import SourceError

        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(side_effect=SourceError("connection error"))
            mock_registry.get.return_value = mock_adapter

            result = await search_async("aspirin")
            assert isinstance(result, CollectionType)
            assert len(result) == 0

    async def test_search_async_query_type_passed(self):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            await search_async("2244", query_type="cid")
            mock_adapter.search.assert_called_once_with("2244", query_type="cid")


class TestSearch:
    def test_search_sync_returns_collection(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=aspirin_collection.compounds)
            mock_registry.get.return_value = mock_adapter

            result = search("aspirin")
            assert isinstance(result, CollectionType)


class TestGetAsync:
    async def test_get_async_returns_compound(self, aspirin_compound: Compound):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_by_id = AsyncMock(return_value=aspirin_compound)
            mock_registry.get.return_value = mock_adapter

            result = await get_async("2244")
            assert result is aspirin_compound

    async def test_get_async_not_found_returns_none(self):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_by_id = AsyncMock(return_value=None)
            mock_registry.get.return_value = mock_adapter

            result = await get_async("99999")
            assert result is None


class TestGet:
    def test_get_sync_returns_compound(self, aspirin_compound: Compound):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_by_id = AsyncMock(return_value=aspirin_compound)
            mock_registry.get.return_value = mock_adapter

            result = get("2244")
            assert result is aspirin_compound


class TestFindSimilarAsync:
    async def test_find_similar_async_returns_collection(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_similarity = AsyncMock(return_value=aspirin_collection.compounds)
            mock_registry.get.return_value = mock_adapter

            result = await find_similar_async("CC(=O)Oc1ccccc1C(=O)O")
            assert isinstance(result, CollectionType)

    async def test_find_similar_async_no_similarity_method(self):
        """If adapter has no get_similarity, should return empty collection."""
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock(spec=[])  # No get_similarity attribute
            mock_registry.get.return_value = mock_adapter

            result = await find_similar_async("CC(=O)Oc1ccccc1C(=O)O")
            assert isinstance(result, CollectionType)
            assert len(result) == 0


class TestFindSimilar:
    def test_find_similar_sync(self, aspirin_collection: CompoundCollection):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.get_similarity = AsyncMock(return_value=aspirin_collection.compounds)
            mock_registry.get.return_value = mock_adapter

            result = find_similar("CC(=O)Oc1ccccc1C(=O)O")
            assert isinstance(result, CollectionType)


class TestCrossReferenceAsync:
    async def test_cross_reference_async_merges_data(self, aspirin_compound: Compound):
        compound_with_extra = Compound(
            smiles=aspirin_compound.smiles,
            inchikey=aspirin_compound.inchikey,
            formula="C9H8O4",
            sources=["pubchem"],
        )
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[compound_with_extra])
            mock_registry.get.return_value = mock_adapter

            result = await chemfuse.cross_reference_async(aspirin_compound)
            assert result is not None

    async def test_cross_reference_no_smiles_skips(self):
        """Compound with no SMILES or inchi should skip cross-referencing."""
        empty_compound = Compound(smiles="", sources=["test"])
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_registry.get.return_value = mock_adapter

            result = await chemfuse.cross_reference_async(empty_compound)
            assert result is not None  # Returns original compound

    async def test_cross_reference_exception_returns_original(self, aspirin_compound: Compound):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(side_effect=Exception("error"))
            mock_registry.get.return_value = mock_adapter

            result = await chemfuse.cross_reference_async(aspirin_compound)
            assert result is not None


class TestCrossReference:
    def test_cross_reference_sync(self, aspirin_compound: Compound):
        with patch("chemfuse.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            result = chemfuse.cross_reference(aspirin_compound)
            assert result is not None


class TestBatchSearchAPI:
    def test_batch_search_delegates_to_core(self, tmp_path):
        query_file = tmp_path / "queries.txt"
        query_file.write_text("aspirin\n")

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[])
            mock_registry.get.return_value = mock_adapter

            collection, errors = chemfuse.batch_search(str(query_file))
            assert isinstance(collection, CollectionType)
            assert isinstance(errors, list)
