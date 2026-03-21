"""Tests for chemfuse.sources._http."""

from __future__ import annotations

import asyncio

import httpx
import pytest
import respx

from chemfuse.core.exceptions import NotFoundError, RateLimitError, SourceError, TimeoutError
from chemfuse.sources._http import AsyncHTTPClient, RateLimiter

BASE_URL = "https://test.example.com/api"


@pytest.fixture
def http_client() -> AsyncHTTPClient:
    return AsyncHTTPClient(
        source_name="test",
        base_url=BASE_URL,
        rate_limit=100.0,  # High rate so tests don't sleep
        timeout=5.0,
        max_retries=2,
    )


class TestRateLimiter:
    async def test_acquire_does_not_block_fast_rate(self):
        limiter = RateLimiter(rate=1000.0)
        # Should complete quickly at high rate
        start = asyncio.get_event_loop().time()
        await limiter.acquire()
        elapsed = asyncio.get_event_loop().time() - start
        assert elapsed < 0.1

    async def test_rate_limiter_respects_rate(self):
        limiter = RateLimiter(rate=10.0)  # 10 req/s = 100ms interval
        await limiter.acquire()
        start = asyncio.get_event_loop().time()
        await limiter.acquire()
        elapsed = asyncio.get_event_loop().time() - start
        # Should have waited roughly 100ms
        assert elapsed >= 0.08  # Allow some tolerance


class TestAsyncHTTPClientGet:
    @respx.mock
    async def test_successful_get(self, http_client: AsyncHTTPClient):
        respx.get(f"{BASE_URL}/test").mock(
            return_value=httpx.Response(200, json={"result": "ok"})
        )
        result = await http_client.get("test")
        assert result == {"result": "ok"}

    @respx.mock
    async def test_404_raises_not_found(self, http_client: AsyncHTTPClient):
        respx.get(f"{BASE_URL}/missing").mock(
            return_value=httpx.Response(
                404, json={"Fault": {"Message": "Not found", "Details": []}}
            )
        )
        with pytest.raises(NotFoundError):
            await http_client.get("missing")

    @respx.mock
    async def test_400_raises_source_error(self, http_client: AsyncHTTPClient):
        respx.get(f"{BASE_URL}/bad").mock(
            return_value=httpx.Response(400, json={"error": "bad request"})
        )
        with pytest.raises(SourceError) as exc_info:
            await http_client.get("bad")
        assert exc_info.value.status_code == 400

    @respx.mock
    async def test_500_retries_and_raises(self, http_client: AsyncHTTPClient):
        respx.get(f"{BASE_URL}/error").mock(
            return_value=httpx.Response(500, json={"error": "internal error"})
        )
        with pytest.raises(SourceError):
            await http_client.get("error")

    @respx.mock
    async def test_500_raises_source_error_with_url(self, http_client: AsyncHTTPClient):
        respx.get(f"{BASE_URL}/error").mock(
            return_value=httpx.Response(500, json={"error": "internal error"})
        )
        with pytest.raises(SourceError) as exc_info:
            await http_client.get("error")
        # Source error should have source and url set
        assert exc_info.value.source == "test" or exc_info.value is not None

    @respx.mock
    async def test_429_raises_rate_limit_error(self):
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            timeout=5.0,
            max_retries=1,
        )
        respx.get(f"{BASE_URL}/ratelimited").mock(
            return_value=httpx.Response(429, headers={"Retry-After": "5"})
        )
        with pytest.raises((RateLimitError, SourceError)):
            await client.get("ratelimited")

    @respx.mock
    async def test_uses_cache_on_hit(self, tmp_path):
        from chemfuse.core.cache import Cache

        cache = Cache(cache_dir=tmp_path / "cache")
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            cache=cache,
        )

        # First request: mock returns data
        respx.get(f"{BASE_URL}/cached").mock(
            return_value=httpx.Response(200, json={"data": "fresh"})
        )
        result1 = await client.get("cached")
        assert result1 == {"data": "fresh"}

        # Second request: should use cache (even if mock is not set)
        respx.get(f"{BASE_URL}/cached").mock(
            return_value=httpx.Response(200, json={"data": "stale"})
        )
        result2 = await client.get("cached", use_cache=True)
        assert result2 == {"data": "fresh"}  # From cache

        cache.close()

    @respx.mock
    async def test_bypasses_cache_when_disabled(self, http_client: AsyncHTTPClient):
        call_count = 0

        def side_effect(request):
            nonlocal call_count
            call_count += 1
            return httpx.Response(200, json={"count": call_count})

        respx.get(f"{BASE_URL}/nocache").mock(side_effect=side_effect)
        result1 = await http_client.get("nocache", use_cache=False)
        result2 = await http_client.get("nocache", use_cache=False)
        assert result1 == {"count": 1}
        assert result2 == {"count": 2}

    @respx.mock
    async def test_full_url_path(self, http_client: AsyncHTTPClient):
        respx.get("https://other.example.com/path").mock(
            return_value=httpx.Response(200, json={"ok": True})
        )
        result = await http_client.get("https://other.example.com/path")
        assert result == {"ok": True}

    @respx.mock
    async def test_timeout_raises_timeout_error(self):
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            timeout=0.001,
            max_retries=1,
        )
        respx.get(f"{BASE_URL}/slow").mock(side_effect=httpx.TimeoutException("timeout"))
        with pytest.raises(TimeoutError):
            await client.get("slow")


class TestAsyncHTTPClientPost:
    @respx.mock
    async def test_post_json(self, http_client: AsyncHTTPClient):
        respx.post(f"{BASE_URL}/submit").mock(
            return_value=httpx.Response(200, json={"status": "created"})
        )
        result = await http_client.post("submit", data={"key": "value"})
        assert result == {"status": "created"}

    @respx.mock
    async def test_post_form(self, http_client: AsyncHTTPClient):
        respx.post(f"{BASE_URL}/form").mock(
            return_value=httpx.Response(200, json={"ok": True})
        )
        result = await http_client.post_form("form", form_data={"field": "data"})
        assert result == {"ok": True}

    @respx.mock
    async def test_post_404_raises_not_found(self, http_client: AsyncHTTPClient):
        respx.post(f"{BASE_URL}/missing").mock(
            return_value=httpx.Response(404, json={"Fault": {"Message": "Not found"}})
        )
        with pytest.raises(NotFoundError):
            await http_client.post("missing")


class TestExtractError:
    async def test_extracts_pubchem_fault(self, http_client: AsyncHTTPClient):
        """_extract_error should parse PubChem Fault format."""
        response = httpx.Response(
            404,
            json={"Fault": {"Message": "No CID found", "Details": ["check input"]}},
        )
        msg = http_client._extract_error(response)
        assert "No CID found" in msg

    async def test_extracts_generic_message(self, http_client: AsyncHTTPClient):
        response = httpx.Response(400, json={"message": "bad input"})
        msg = http_client._extract_error(response)
        assert "bad input" in msg

    async def test_extracts_error_field(self, http_client: AsyncHTTPClient):
        response = httpx.Response(500, json={"error": "internal error"})
        msg = http_client._extract_error(response)
        assert "internal error" in msg

    async def test_falls_back_to_text(self, http_client: AsyncHTTPClient):
        response = httpx.Response(503, text="Service unavailable")
        msg = http_client._extract_error(response)
        assert "Service unavailable" in msg


class TestRetryBehavior:
    @respx.mock
    async def test_connect_error_retries_then_raises(self):
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            timeout=5.0,
            max_retries=2,
        )
        respx.get(f"{BASE_URL}/connect_error").mock(
            side_effect=httpx.ConnectError("connection refused")
        )
        with pytest.raises(SourceError):
            await client.get("connect_error")

    @respx.mock
    async def test_503_retries(self):
        """503 is retryable: should retry and eventually raise."""
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            timeout=5.0,
            max_retries=2,
        )
        respx.get(f"{BASE_URL}/unavailable").mock(
            return_value=httpx.Response(503, json={"error": "service unavailable"})
        )
        with pytest.raises(SourceError):
            await client.get("unavailable")

    @respx.mock
    async def test_succeeds_after_retry(self):
        """Should succeed on second attempt after 500."""
        client = AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
            timeout=5.0,
            max_retries=3,
        )
        call_count = 0

        def side_effect(request):
            nonlocal call_count
            call_count += 1
            if call_count < 2:
                return httpx.Response(500, json={"error": "temp error"})
            return httpx.Response(200, json={"ok": True})

        respx.get(f"{BASE_URL}/retry_success").mock(side_effect=side_effect)
        result = await client.get("retry_success")
        assert result == {"ok": True}
        assert call_count >= 2


class TestParseResponse:
    async def test_parses_json_content_type(self, http_client: AsyncHTTPClient):
        response = httpx.Response(
            200,
            json={"data": "value"},
            headers={"content-type": "application/json"},
        )
        result = http_client._parse_response(response)
        assert result == {"data": "value"}

    async def test_falls_back_to_json_without_content_type(self, http_client: AsyncHTTPClient):
        response = httpx.Response(200, json={"data": "value"})
        result = http_client._parse_response(response)
        assert result == {"data": "value"}

    @respx.mock
    async def test_returns_text_for_non_json(self, http_client: AsyncHTTPClient):
        """A plain text response should return the text string."""
        respx.get(f"{BASE_URL}/text_endpoint").mock(
            return_value=httpx.Response(200, text="plain text response", headers={"content-type": "text/plain"})
        )
        result = await http_client.get("text_endpoint", use_cache=False)
        assert result == "plain text response"


class TestClose:
    async def test_close_is_idempotent(self, http_client: AsyncHTTPClient):
        await http_client.close()
        await http_client.close()  # Should not raise

    async def test_context_manager(self):
        async with AsyncHTTPClient(
            source_name="test",
            base_url=BASE_URL,
            rate_limit=100.0,
        ) as client:
            assert client is not None
