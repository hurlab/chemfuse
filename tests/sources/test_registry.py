"""Tests for chemfuse.sources SourceRegistry."""

from __future__ import annotations

import pytest

from chemfuse.models.compound import Compound
from chemfuse.sources import SourceRegistry
from chemfuse.sources._base import SourceAdapter


class MockAdapter(SourceAdapter):
    name = "mock"
    base_url = "https://mock.example.com"
    rate_limit = 10.0

    async def search(self, query: str, query_type: str = "name") -> list[Compound]:
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        return None

    async def get_properties(self, identifier: str) -> dict:
        return {}

    def is_available(self) -> bool:
        return True


class AnotherAdapter(SourceAdapter):
    name = "another"
    base_url = "https://another.example.com"
    rate_limit = 5.0

    async def search(self, query: str, query_type: str = "name") -> list[Compound]:
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        return None

    async def get_properties(self, identifier: str) -> dict:
        return {}

    def is_available(self) -> bool:
        return True


@pytest.fixture
def registry() -> SourceRegistry:
    """Fresh SourceRegistry for each test."""
    return SourceRegistry()


class TestSourceRegistry:
    def test_register_and_list(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        assert "mock" in registry.list()

    def test_list_sorted(self, registry: SourceRegistry):
        registry.register("zebra", MockAdapter)
        registry.register("apple", AnotherAdapter)
        names = registry.list()
        assert names == sorted(names)

    def test_get_instantiates_adapter(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        adapter = registry.get("mock")
        assert isinstance(adapter, MockAdapter)

    def test_get_returns_same_instance(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        adapter1 = registry.get("mock")
        adapter2 = registry.get("mock")
        assert adapter1 is adapter2

    def test_get_unknown_raises_key_error(self, registry: SourceRegistry):
        with pytest.raises(KeyError) as exc_info:
            registry.get("nonexistent")
        assert "nonexistent" in str(exc_info.value)

    def test_key_error_message_shows_available(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        with pytest.raises(KeyError) as exc_info:
            registry.get("nonexistent")
        assert "mock" in str(exc_info.value)

    def test_is_registered_true(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        assert registry.is_registered("mock") is True

    def test_is_registered_false(self, registry: SourceRegistry):
        assert registry.is_registered("nonexistent") is False

    def test_reregister_clears_cached_instance(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        instance1 = registry.get("mock")
        registry.register("mock", MockAdapter)  # Re-register same class
        instance2 = registry.get("mock")
        assert instance1 is not instance2

    def test_empty_registry_list(self, registry: SourceRegistry):
        assert registry.list() == []

    def test_multiple_adapters(self, registry: SourceRegistry):
        registry.register("mock", MockAdapter)
        registry.register("another", AnotherAdapter)
        assert len(registry.list()) == 2
        assert registry.is_registered("mock")
        assert registry.is_registered("another")


class TestGlobalRegistry:
    def test_global_registry_has_pubchem(self):
        from chemfuse.sources import registry
        assert registry.is_registered("pubchem")

    def test_global_registry_returns_pubchem_adapter(self):
        from chemfuse.sources import registry
        from chemfuse.sources.pubchem import PubChemAdapter

        adapter = registry.get("pubchem")
        assert isinstance(adapter, PubChemAdapter)
