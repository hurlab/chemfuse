"""Tests for chemfuse.core.config."""

import os

from chemfuse.core.config import Settings, get_settings, reset_settings


class TestSettings:
    def setup_method(self):
        reset_settings()

    def teardown_method(self):
        reset_settings()

    def test_default_values(self):
        settings = Settings()
        assert settings.cache_enabled is True
        assert settings.cache_ttl == 604800
        assert settings.request_timeout == 30.0
        assert settings.max_retries == 3
        assert "pubchem" in settings.default_sources

    def test_cache_dir_default(self):
        settings = Settings()
        assert ".chemfuse" in str(settings.cache_dir)

    def test_rate_limits_default(self):
        settings = Settings()
        assert "pubchem" in settings.rate_limits
        assert settings.rate_limits["pubchem"] > 0

    def test_env_prefix(self):
        """Settings should use CHEMFUSE_ prefix for env vars."""
        os.environ["CHEMFUSE_REQUEST_TIMEOUT"] = "60.0"
        try:
            settings = Settings()
            assert settings.request_timeout == 60.0
        finally:
            del os.environ["CHEMFUSE_REQUEST_TIMEOUT"]

    def test_cache_disabled_via_env(self):
        os.environ["CHEMFUSE_CACHE_ENABLED"] = "false"
        try:
            settings = Settings()
            assert settings.cache_enabled is False
        finally:
            del os.environ["CHEMFUSE_CACHE_ENABLED"]

    def test_max_retries_configurable(self):
        os.environ["CHEMFUSE_MAX_RETRIES"] = "5"
        try:
            settings = Settings()
            assert settings.max_retries == 5
        finally:
            del os.environ["CHEMFUSE_MAX_RETRIES"]


class TestGetSettings:
    def setup_method(self):
        reset_settings()

    def teardown_method(self):
        reset_settings()

    def test_returns_settings_instance(self):
        s = get_settings()
        assert isinstance(s, Settings)

    def test_singleton_behavior(self):
        s1 = get_settings()
        s2 = get_settings()
        assert s1 is s2

    def test_reset_creates_new_instance(self):
        s1 = get_settings()
        reset_settings()
        s2 = get_settings()
        assert s1 is not s2

    def test_settings_are_correct_after_reset(self):
        s = get_settings()
        assert s.cache_enabled is True
        reset_settings()
        s2 = get_settings()
        assert s2.cache_enabled is True
