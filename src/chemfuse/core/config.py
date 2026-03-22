"""ChemFuse configuration management using pydantic-settings."""

from __future__ import annotations

import threading
from pathlib import Path
from typing import Any

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """ChemFuse configuration settings.

    Settings can be overridden via environment variables prefixed with CHEMFUSE_.
    For example: CHEMFUSE_CACHE_DIR=/tmp/cache
    """

    model_config = SettingsConfigDict(
        env_prefix="CHEMFUSE_",
        env_file=None,
        case_sensitive=False,
    )

    # Cache settings
    cache_dir: Path = Field(
        default_factory=lambda: Path.home() / ".chemfuse" / "cache",
        description="Directory for SQLite cache database",
    )
    cache_enabled: bool = Field(default=True, description="Enable response caching")
    cache_ttl: int = Field(
        default=604800,  # 7 days
        description="Cache time-to-live in seconds",
    )
    cache_max_entries: int = Field(
        default=10000,
        description="Maximum number of entries in the cache (0 = unlimited). When exceeded, the oldest 10% of entries (by last_accessed) are evicted.",
    )

    # Source settings
    default_sources: list[str] = Field(
        default_factory=lambda: ["pubchem"],
        description="Default data sources to use",
    )

    # Timeout settings (seconds)
    request_timeout: float = Field(default=30.0, description="HTTP request timeout")
    list_key_timeout: float = Field(
        default=120.0, description="Timeout for ListKey polling"
    )

    # Rate limits per source (requests per second)
    rate_limits: dict[str, float] = Field(
        default_factory=lambda: {
            "pubchem": 5.0,
            "chembl": 1.0,
            "unichem": 5.0,
            "bindingdb": 1.0,
            "opentargets": 5.0,
            "surechembl": 2.0,
        },
        description="Per-source rate limits in requests/second",
    )

    # Retry settings
    max_retries: int = Field(default=3, description="Maximum HTTP retry attempts")

    @field_validator("cache_dir", mode="before")
    @classmethod
    def expand_cache_dir(cls, v: Any) -> Path:
        """Expand ~ in cache directory path."""
        return Path(str(v)).expanduser()

    def get_rate_limit(self, source: str) -> float:
        """Get rate limit for a specific source.

        Args:
            source: Source name.

        Returns:
            Rate limit in requests per second.
        """
        return self.rate_limits.get(source, 1.0)


# Module-level singleton with thread safety
_settings: Settings | None = None
_lock = threading.Lock()


def get_settings() -> Settings:
    """Get or create the global Settings instance.

    Uses double-checked locking for thread-safe singleton initialization.

    Returns:
        The global Settings instance.
    """
    global _settings
    if _settings is None:
        with _lock:
            if _settings is None:
                _settings = Settings()
    return _settings


def reset_settings() -> None:
    """Reset the global Settings instance (for testing)."""
    global _settings
    _settings = None
