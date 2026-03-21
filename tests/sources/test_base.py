"""Tests for chemfuse.sources._base."""

from __future__ import annotations

import pytest

from chemfuse.models.compound import Compound
from chemfuse.sources._base import SourceAdapter


class ConcreteAdapter(SourceAdapter):
    """Minimal concrete implementation for testing."""

    name = "concrete"
    base_url = "https://example.com"
    rate_limit = 1.0

    async def search(self, query: str, query_type: str = "name") -> list[Compound]:
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        return None

    async def get_properties(self, identifier: str) -> dict:
        return {}

    def is_available(self) -> bool:
        return True


class IncompleteAdapter(SourceAdapter):
    """Missing required methods."""
    name = "incomplete"


class TestSourceAdapterABC:
    def test_cannot_instantiate_abstract(self):
        with pytest.raises(TypeError):
            SourceAdapter()  # type: ignore[abstract]

    def test_cannot_instantiate_incomplete(self):
        with pytest.raises(TypeError):
            IncompleteAdapter()  # type: ignore[abstract]

    def test_concrete_adapter_instantiates(self):
        adapter = ConcreteAdapter()
        assert adapter is not None

    def test_name_attribute(self):
        adapter = ConcreteAdapter()
        assert adapter.name == "concrete"

    def test_base_url_attribute(self):
        adapter = ConcreteAdapter()
        assert adapter.base_url == "https://example.com"

    def test_rate_limit_attribute(self):
        adapter = ConcreteAdapter()
        assert adapter.rate_limit == 1.0

    async def test_search_returns_list(self):
        adapter = ConcreteAdapter()
        result = await adapter.search("aspirin")
        assert isinstance(result, list)

    async def test_get_by_id_returns_none(self):
        adapter = ConcreteAdapter()
        result = await adapter.get_by_id("12345")
        assert result is None

    async def test_get_properties_returns_dict(self):
        adapter = ConcreteAdapter()
        result = await adapter.get_properties("12345")
        assert isinstance(result, dict)

    def test_is_available_returns_bool(self):
        adapter = ConcreteAdapter()
        assert isinstance(adapter.is_available(), bool)
