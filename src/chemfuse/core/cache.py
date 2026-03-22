"""SQLite-based response cache with TTL support."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import time
from pathlib import Path
from typing import Any

DEFAULT_CACHE_DIR = Path.home() / ".chemfuse" / "cache"
DEFAULT_TTL_SECONDS = 604800  # 7 days


class Cache:
    """SQLite-based cache for HTTP API responses.

    Stores JSON-serializable responses keyed by SHA-256 hash of URL and
    sorted parameters, with configurable TTL expiration and hit counting.
    """

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        ttl: int = DEFAULT_TTL_SECONDS,
        enabled: bool = True,
    ) -> None:
        """Initialize the cache.

        Args:
            cache_dir: Directory to store the SQLite database. Defaults to ~/.chemfuse/cache.
            ttl: Time-to-live in seconds for cached entries. Defaults to 7 days.
            enabled: Whether caching is enabled.
        """
        self.enabled = enabled
        self.ttl = ttl

        if cache_dir is None:
            cache_dir = DEFAULT_CACHE_DIR
        self.cache_dir = Path(cache_dir).expanduser()
        self._db_path = self.cache_dir / "chemfuse_cache.db"
        self._conn: sqlite3.Connection | None = None

        if self.enabled:
            self._init_db()

    def _init_db(self) -> None:
        """Initialize the SQLite database and create tables if needed."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self._db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS cache (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL,
                created_at REAL NOT NULL,
                expires_at REAL NOT NULL,
                url TEXT,
                hit_count INTEGER DEFAULT 0
            )
        """)
        self._conn.execute("""
            CREATE INDEX IF NOT EXISTS idx_expires_at ON cache (expires_at)
        """)
        self._conn.commit()

    @staticmethod
    def _make_key(url: str, params: dict[str, Any] | None = None) -> str:
        """Generate a cache key from URL and parameters.

        Args:
            url: The request URL.
            params: Optional request parameters.

        Returns:
            A SHA-256 hash string as the cache key.
        """
        key_data = url
        if params:
            key_data += json.dumps(params, sort_keys=True)
        return hashlib.sha256(key_data.encode()).hexdigest()

    def get(self, url: str, params: dict[str, Any] | None = None) -> Any | None:
        """Retrieve a cached response.

        Args:
            url: The request URL.
            params: Optional request parameters.

        Returns:
            The cached value if found and not expired, otherwise None.
        """
        if not self.enabled or self._conn is None:
            return None

        key = self._make_key(url, params)
        now = time.time()

        cursor = self._conn.execute(
            "SELECT value FROM cache WHERE key = ? AND expires_at > ?",
            (key, now),
        )
        row = cursor.fetchone()
        if row is not None:
            self._conn.execute(
                "UPDATE cache SET hit_count = hit_count + 1 WHERE key = ?",
                (key,),
            )
            self._conn.commit()
            return json.loads(row[0])
        return None

    def set(
        self,
        url: str,
        value: Any,
        params: dict[str, Any] | None = None,
        ttl: int | None = None,
    ) -> None:
        """Store a response in the cache.

        Args:
            url: The request URL.
            value: The response value (must be JSON-serializable).
            params: Optional request parameters.
            ttl: Optional TTL override in seconds.
        """
        if not self.enabled or self._conn is None:
            return

        key = self._make_key(url, params)
        now = time.time()
        expires_at = now + (ttl if ttl is not None else self.ttl)

        self._conn.execute(
            """INSERT OR REPLACE INTO cache (key, value, created_at, expires_at, url, hit_count)
               VALUES (?, ?, ?, ?, ?, 0)""",
            (key, json.dumps(value), now, expires_at, url),
        )
        self._conn.commit()

    def invalidate(self, url: str, params: dict[str, Any] | None = None) -> bool:
        """Delete a specific cached entry.

        Args:
            url: The request URL.
            params: Optional request parameters.

        Returns:
            True if an entry was deleted.
        """
        if not self.enabled or self._conn is None:
            return False

        key = self._make_key(url, params)
        cursor = self._conn.execute("DELETE FROM cache WHERE key = ?", (key,))
        self._conn.commit()
        return cursor.rowcount > 0

    def clear(self) -> int:
        """Remove all cached entries.

        Returns:
            Number of entries removed.
        """
        if not self.enabled or self._conn is None:
            return 0

        cursor = self._conn.execute("DELETE FROM cache")
        self._conn.commit()
        return cursor.rowcount

    def clear_expired(self) -> int:
        """Remove all expired entries.

        Returns:
            Number of entries removed.
        """
        if not self.enabled or self._conn is None:
            return 0

        now = time.time()
        cursor = self._conn.execute("DELETE FROM cache WHERE expires_at <= ?", (now,))
        self._conn.commit()
        return cursor.rowcount

    def stats(self) -> dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary with total_entries, expired_entries, total_hits, db_size_bytes.
        """
        if not self.enabled or self._conn is None:
            return {"total_entries": 0, "expired_entries": 0, "total_hits": 0, "db_size_bytes": 0}

        now = time.time()
        total = self._conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
        expired = self._conn.execute(
            "SELECT COUNT(*) FROM cache WHERE expires_at <= ?", (now,)
        ).fetchone()[0]
        hits = self._conn.execute(
            "SELECT COALESCE(SUM(hit_count), 0) FROM cache"
        ).fetchone()[0]
        db_size = self._db_path.stat().st_size if self._db_path.exists() else 0

        return {
            "total_entries": total,
            "expired_entries": expired,
            "total_hits": hits,
            "db_size_bytes": db_size,
        }

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> Cache:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()
