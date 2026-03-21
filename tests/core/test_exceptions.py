"""Tests for chemfuse.core.exceptions."""

import pytest

from chemfuse.core.exceptions import (
    ChemFuseError,
    NotFoundError,
    RateLimitError,
    SourceError,
    TimeoutError,
    ValidationError,
)


class TestChemFuseError:
    def test_base_error_is_exception(self):
        err = ChemFuseError("test message")
        assert isinstance(err, Exception)
        assert str(err) == "test message"

    def test_base_error_subclass(self):
        assert issubclass(SourceError, ChemFuseError)
        assert issubclass(NotFoundError, ChemFuseError)
        assert issubclass(RateLimitError, ChemFuseError)
        assert issubclass(ValidationError, ChemFuseError)
        assert issubclass(TimeoutError, ChemFuseError)


class TestSourceError:
    def test_basic_construction(self):
        err = SourceError("pubchem error")
        assert "pubchem error" in str(err)

    def test_with_source_and_url(self):
        err = SourceError("request failed", source="pubchem", url="https://example.com", status_code=500)
        assert err.source == "pubchem"
        assert err.url == "https://example.com"
        assert err.status_code == 500

    def test_defaults(self):
        err = SourceError("msg")
        assert err.source is None
        assert err.url is None
        assert err.status_code is None

    def test_raise_and_catch(self):
        with pytest.raises(SourceError) as exc_info:
            raise SourceError("test", source="pubchem", status_code=500)
        assert exc_info.value.source == "pubchem"
        assert exc_info.value.status_code == 500


class TestNotFoundError:
    def test_basic(self):
        err = NotFoundError("compound not found")
        assert "compound not found" in str(err)

    def test_with_identifier(self):
        err = NotFoundError("not found", identifier="aspirin")
        assert err.identifier == "aspirin"

    def test_is_chemfuse_error(self):
        with pytest.raises(ChemFuseError):
            raise NotFoundError("test")

    def test_default_identifier(self):
        err = NotFoundError("msg")
        assert err.identifier is None


class TestRateLimitError:
    def test_basic(self):
        err = RateLimitError("rate limited")
        assert "rate limited" in str(err)

    def test_with_retry_after(self):
        err = RateLimitError("retry after 5s", retry_after=5.0)
        assert err.retry_after == 5.0

    def test_default_retry_after(self):
        err = RateLimitError("msg")
        assert err.retry_after is None


class TestValidationError:
    def test_basic(self):
        err = ValidationError("invalid input")
        assert "invalid input" in str(err)

    def test_is_chemfuse_error(self):
        assert isinstance(ValidationError("msg"), ChemFuseError)


class TestTimeoutError:
    def test_basic(self):
        err = TimeoutError("request timed out")
        assert "request timed out" in str(err)

    def test_is_chemfuse_error(self):
        assert isinstance(TimeoutError("msg"), ChemFuseError)

    def test_raise_and_catch_as_base(self):
        with pytest.raises(ChemFuseError):
            raise TimeoutError("timeout")
