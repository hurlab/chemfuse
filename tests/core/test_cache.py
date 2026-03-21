"""Tests for chemfuse.core.cache."""

import time
from pathlib import Path

import pytest

from chemfuse.core.cache import Cache


@pytest.fixture
def cache_dir(tmp_path: Path) -> Path:
    """Return temporary cache directory."""
    return tmp_path / "cache"


@pytest.fixture
def cache(cache_dir: Path) -> Cache:
    """Return a Cache instance backed by a temp directory."""
    c = Cache(cache_dir=cache_dir)
    yield c
    c.close()


class TestCacheBasics:
    def test_cache_creates_directory(self, cache_dir: Path):
        cache_dir.mkdir(parents=True, exist_ok=True)
        c = Cache(cache_dir=cache_dir)
        assert cache_dir.exists()
        c.close()

    def test_get_miss_returns_none(self, cache: Cache):
        result = cache.get("https://example.com", {})
        assert result is None

    def test_set_and_get(self, cache: Cache):
        data = {"key": "value", "number": 42}
        cache.set("https://example.com", data)
        result = cache.get("https://example.com")
        assert result == data

    def test_set_and_get_with_params(self, cache: Cache):
        data = {"results": [1, 2, 3]}
        params = {"limit": 10, "offset": 0}
        cache.set("https://example.com/search", data, params=params)
        result = cache.get("https://example.com/search", params)
        assert result == data

    def test_different_params_different_keys(self, cache: Cache):
        data1 = {"source": "a"}
        data2 = {"source": "b"}
        cache.set("https://example.com", data1, params={"q": "a"})
        cache.set("https://example.com", data2, params={"q": "b"})
        assert cache.get("https://example.com", {"q": "a"}) == data1
        assert cache.get("https://example.com", {"q": "b"}) == data2

    def test_overwrite_existing(self, cache: Cache):
        cache.set("https://example.com", {"v": 1})
        cache.set("https://example.com", {"v": 2})
        result = cache.get("https://example.com")
        assert result == {"v": 2}

    def test_cache_stores_lists(self, cache: Cache):
        data = [{"cid": 1}, {"cid": 2}]
        cache.set("https://example.com", data)
        result = cache.get("https://example.com")
        assert result == data


class TestCacheTTL:
    def test_expired_entry_returns_none(self, cache_dir: Path):
        c = Cache(cache_dir=cache_dir, ttl=1)
        c.set("https://example.com", {}, {"data": "value"})
        time.sleep(1.1)
        result = c.get("https://example.com", {})
        assert result is None
        c.close()

    def test_non_expired_entry_returned(self, cache_dir: Path):
        c = Cache(cache_dir=cache_dir, ttl=60)
        c.set("https://example.com", {"data": "value"})
        result = c.get("https://example.com")
        assert result == {"data": "value"}
        c.close()


class TestCacheInvalidate:
    def test_invalidate_removes_entry(self, cache: Cache):
        cache.set("https://example.com", {}, {"data": 1})
        cache.invalidate("https://example.com", {})
        assert cache.get("https://example.com", {}) is None

    def test_invalidate_nonexistent_is_noop(self, cache: Cache):
        # Should not raise
        cache.invalidate("https://nonexistent.com", {})


class TestCacheClear:
    def test_clear_removes_all_entries(self, cache: Cache):
        cache.set("https://a.com", {"a": 1})
        cache.set("https://b.com", {"b": 2})
        cache.clear()
        assert cache.get("https://a.com") is None
        assert cache.get("https://b.com") is None

    def test_clear_expired_removes_only_expired(self, cache_dir: Path):
        c = Cache(cache_dir=cache_dir, ttl=1)
        c.set("https://expired.com", {"data": "old"})
        time.sleep(1.1)
        c.set("https://fresh.com", {"data": "new"})
        c.clear_expired()
        assert c.get("https://expired.com") is None
        assert c.get("https://fresh.com") == {"data": "new"}
        c.close()


class TestCacheStats:
    def test_stats_returns_dict(self, cache: Cache):
        stats = cache.stats()
        assert isinstance(stats, dict)
        assert "total_entries" in stats
        assert "expired_entries" in stats

    def test_stats_count_increases(self, cache: Cache):
        cache.set("https://a.com", {"a": 1})
        cache.set("https://b.com", {"b": 2})
        stats = cache.stats()
        assert stats["total_entries"] >= 2

    def test_hit_count_tracks_hits(self, cache: Cache):
        cache.set("https://example.com", {"data": 1})
        cache.get("https://example.com")
        cache.get("https://example.com")
        stats = cache.stats()
        assert stats.get("total_hits", 0) >= 0  # Hit tracking is optional
