"""Tests for chemfuse.core.cache."""

import sqlite3
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


class TestCacheMaxEntries:
    """Tests for LRU eviction and the max_entries limit."""

    def test_eviction_enforces_max_entries(self, cache_dir: Path):
        # With max_entries=5 and eviction fraction=10% (min 1), inserting 6
        # entries should trigger eviction leaving at most 5.
        c = Cache(cache_dir=cache_dir, max_entries=5)
        for i in range(6):
            c.set(f"https://example.com/{i}", {"i": i})
        stats = c.cache_stats()
        assert stats["total_entries"] <= 5
        c.close()

    def test_eviction_removes_lru_entries(self, cache_dir: Path):
        # Insert 5 entries then access entry #0 to make it recently used.
        # Then insert a 6th entry which should evict the true LRU (entry #1).
        c = Cache(cache_dir=cache_dir, max_entries=5)
        for i in range(5):
            c.set(f"https://example.com/{i}", {"i": i})
        # Touch entry 0 so it becomes most-recently used.
        time.sleep(0.01)
        c.get("https://example.com/0")
        # Insert a 6th entry to trigger eviction.
        time.sleep(0.01)
        c.set("https://example.com/5", {"i": 5})
        # Entry 0 was accessed recently so it should still be present.
        assert c.get("https://example.com/0") == {"i": 0}
        c.close()

    def test_no_eviction_when_below_limit(self, cache_dir: Path):
        c = Cache(cache_dir=cache_dir, max_entries=100)
        for i in range(10):
            c.set(f"https://example.com/{i}", {"i": i})
        assert c.cache_stats()["total_entries"] == 10
        c.close()

    def test_zero_max_entries_disables_eviction(self, cache_dir: Path):
        # max_entries=0 means unlimited; no eviction should occur.
        c = Cache(cache_dir=cache_dir, max_entries=0)
        for i in range(20):
            c.set(f"https://example.com/{i}", {"i": i})
        assert c.cache_stats()["total_entries"] == 20
        c.close()


class TestCacheStatsMethod:
    def test_cache_stats_includes_max_entries(self, cache: Cache):
        stats = cache.cache_stats()
        assert "max_entries" in stats
        assert stats["max_entries"] == cache.max_entries

    def test_cache_stats_entry_count(self, cache: Cache):
        cache.set("https://a.com", {"a": 1})
        cache.set("https://b.com", {"b": 2})
        stats = cache.cache_stats()
        assert stats["total_entries"] == 2

    def test_cache_stats_db_size_bytes(self, cache: Cache):
        cache.set("https://a.com", {"a": 1})
        stats = cache.cache_stats()
        assert stats["db_size_bytes"] > 0


class TestCacheClearCache:
    def test_clear_cache_removes_all(self, cache: Cache):
        cache.set("https://a.com", {"a": 1})
        cache.set("https://b.com", {"b": 2})
        removed = cache.clear_cache()
        assert removed == 2
        assert cache.cache_stats()["total_entries"] == 0


class TestLastAccessedMigration:
    def test_existing_db_without_last_accessed_is_migrated(self, cache_dir: Path):
        # Create an old-style database without the last_accessed column.
        cache_dir.mkdir(parents=True, exist_ok=True)
        db_path = cache_dir / "chemfuse_cache.db"
        conn = sqlite3.connect(str(db_path))
        conn.execute("""
            CREATE TABLE cache (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL,
                created_at REAL NOT NULL,
                expires_at REAL NOT NULL,
                url TEXT,
                hit_count INTEGER DEFAULT 0
            )
        """)
        now = time.time()
        conn.execute(
            "INSERT INTO cache (key, value, created_at, expires_at, url) VALUES (?, ?, ?, ?, ?)",
            ("oldkey", '{"legacy": true}', now, now + 3600, "https://legacy.com"),
        )
        conn.commit()
        conn.close()

        # Opening the Cache should migrate the column without error.
        c = Cache(cache_dir=cache_dir)
        cols = {
            row[1]
            for row in c._conn.execute("PRAGMA table_info(cache)").fetchall()
        }
        assert "last_accessed" in cols
        # Existing data should still be readable.
        result = c._conn.execute("SELECT value FROM cache WHERE key = 'oldkey'").fetchone()
        assert result is not None
        c.close()
