"""Shared async HTTP client with per-source rate limiting and retry logic."""

from __future__ import annotations

import asyncio
import time
from typing import TYPE_CHECKING, Any

import httpx

from chemfuse.core.exceptions import NotFoundError, RateLimitError, SourceError, TimeoutError

if TYPE_CHECKING:
    from chemfuse.core.cache import Cache


class RateLimiter:
    """Token bucket rate limiter for a single source."""

    def __init__(self, rate: float) -> None:
        """Initialize with requests per second rate.

        Args:
            rate: Maximum requests per second.

        Raises:
            ValueError: If rate is not positive.
        """
        if rate <= 0:
            raise ValueError(f"Rate must be positive, got {rate}")
        self.rate = rate
        self._min_interval = 1.0 / rate
        self._last_request_time: float = 0.0
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        """Acquire permission to make a request, sleeping if necessary."""
        async with self._lock:
            now = time.monotonic()
            elapsed = now - self._last_request_time
            if elapsed < self._min_interval:
                await asyncio.sleep(self._min_interval - elapsed)
            self._last_request_time = time.monotonic()


class AsyncHTTPClient:
    """Shared async HTTP client with rate limiting, caching, and retry.

    Supports per-source rate limiting via RateLimiter instances.
    Implements exponential backoff retries for transient errors.
    """

    RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}

    def __init__(
        self,
        source_name: str,
        base_url: str,
        rate_limit: float = 1.0,
        timeout: float = 30.0,
        max_retries: int = 3,
        cache: Cache | None = None,
        headers: dict[str, str] | None = None,
    ) -> None:
        """Initialize the HTTP client.

        Args:
            source_name: Name of the data source (for error messages).
            base_url: Base URL for the API.
            rate_limit: Maximum requests per second.
            timeout: Request timeout in seconds.
            max_retries: Maximum number of retry attempts.
            cache: Optional Cache instance.
            headers: Optional default headers.
        """
        self.source_name = source_name
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.max_retries = max_retries
        self.cache = cache
        self._rate_limiter = RateLimiter(rate_limit)
        self._client: httpx.AsyncClient | None = None
        self._default_headers = {
            "Accept": "application/json",
            "User-Agent": "ChemFuse/0.1.0 (Python; httpx)",
            **(headers or {}),
        }

    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create the httpx async client."""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                timeout=httpx.Timeout(self.timeout),
                headers=self._default_headers,
                follow_redirects=True,
            )
        return self._client

    async def get(
        self,
        path: str,
        params: dict[str, Any] | None = None,
        use_cache: bool = True,
    ) -> Any:
        """Make a GET request.

        Args:
            path: URL path relative to base_url, or full URL.
            params: Query parameters.
            use_cache: Whether to use the cache.

        Returns:
            Parsed JSON response.

        Raises:
            SourceError: On API errors.
            NotFoundError: When resource is not found (404).
            RateLimitError: When rate limit is exhausted after retries.
            TimeoutError: When request times out after retries.
        """
        url = path if path.startswith("http") else f"{self.base_url}/{path.lstrip('/')}"

        # Check cache
        if use_cache and self.cache:
            cached = self.cache.get(url, params)
            if cached is not None:
                return cached

        last_exc: Exception | None = None
        for attempt in range(self.max_retries):
            await self._rate_limiter.acquire()
            try:
                client = await self._get_client()
                response = await client.get(url, params=params)

                if response.status_code == 200:
                    data = self._parse_response(response)
                    if use_cache and self.cache:
                        self.cache.set(url, data, params)
                    return data

                elif response.status_code == 404:
                    error_msg = self._extract_error(response)
                    raise NotFoundError(
                        f"Not found: {error_msg}",
                    )

                elif response.status_code == 400:
                    error_msg = self._extract_error(response)
                    raise SourceError(
                        f"Bad request: {error_msg}",
                        source=self.source_name,
                        url=url,
                        status_code=400,
                    )

                elif response.status_code in self.RETRYABLE_STATUS_CODES:
                    wait_time = (2.0 ** attempt)
                    await asyncio.sleep(wait_time)
                    if response.status_code == 429:
                        retry_after_str = response.headers.get("Retry-After")
                        retry_after = float(retry_after_str) if retry_after_str else None
                        last_exc = RateLimitError(
                            f"Rate limit exceeded for {self.source_name}",
                            retry_after=retry_after,
                        )
                    else:
                        error_msg = self._extract_error(response)
                        last_exc = SourceError(
                            f"Server error: {error_msg}",
                            source=self.source_name,
                            url=url,
                            status_code=response.status_code,
                        )
                    continue

                else:
                    error_msg = self._extract_error(response)
                    raise SourceError(
                        f"HTTP {response.status_code}: {error_msg}",
                        source=self.source_name,
                        url=url,
                        status_code=response.status_code,
                    )

            except (NotFoundError, SourceError):
                raise
            except httpx.TimeoutException as exc:
                last_exc = TimeoutError(
                    f"Request timed out for {self.source_name}: {url}"
                )
                if attempt < self.max_retries - 1:
                    await asyncio.sleep(1.0 * (attempt + 1))
                    continue
                raise TimeoutError(
                    f"Request timed out after {self.max_retries} attempts: {url}"
                ) from exc
            except httpx.ConnectError as exc:
                last_exc = SourceError(
                    f"Connection error for {self.source_name}",
                    source=self.source_name,
                    url=url,
                )
                if attempt < self.max_retries - 1:
                    await asyncio.sleep(1.0 * (attempt + 1))
                    continue
                raise SourceError(
                    f"Connection failed after {self.max_retries} attempts",
                    source=self.source_name,
                    url=url,
                ) from exc

        if last_exc:
            raise last_exc
        raise SourceError(
            "Request failed with no error details",
            source=self.source_name,
            url=url,
        )

    async def post(
        self,
        path: str,
        data: dict[str, Any] | None = None,
        params: dict[str, Any] | None = None,
    ) -> Any:
        """Make a POST request.

        Args:
            path: URL path relative to base_url, or full URL.
            data: POST body (sent as JSON).
            params: Query parameters.

        Returns:
            Parsed JSON response.
        """
        url = path if path.startswith("http") else f"{self.base_url}/{path.lstrip('/')}"
        await self._rate_limiter.acquire()

        try:
            client = await self._get_client()
            response = await client.post(url, json=data, params=params)

            if response.status_code == 200:
                return self._parse_response(response)
            elif response.status_code == 404:
                error_msg = self._extract_error(response)
                raise NotFoundError(f"Not found: {error_msg}")
            else:
                error_msg = self._extract_error(response)
                raise SourceError(
                    f"HTTP {response.status_code}: {error_msg}",
                    source=self.source_name,
                    url=url,
                    status_code=response.status_code,
                )
        except (NotFoundError, SourceError):
            raise
        except httpx.TimeoutException as exc:
            raise TimeoutError(
                f"POST request timed out for {self.source_name}"
            ) from exc

    async def post_form(
        self,
        path: str,
        form_data: dict[str, Any],
        params: dict[str, Any] | None = None,
    ) -> Any:
        """Make a POST request with form-encoded data.

        Args:
            path: URL path relative to base_url, or full URL.
            form_data: Form data dictionary.
            params: Query parameters.

        Returns:
            Parsed JSON response.
        """
        url = path if path.startswith("http") else f"{self.base_url}/{path.lstrip('/')}"
        await self._rate_limiter.acquire()

        try:
            client = await self._get_client()
            response = await client.post(url, data=form_data, params=params)

            if response.status_code == 200:
                return self._parse_response(response)
            elif response.status_code == 404:
                error_msg = self._extract_error(response)
                raise NotFoundError(f"Not found: {error_msg}")
            else:
                error_msg = self._extract_error(response)
                raise SourceError(
                    f"HTTP {response.status_code}: {error_msg}",
                    source=self.source_name,
                    url=url,
                    status_code=response.status_code,
                )
        except (NotFoundError, SourceError):
            raise
        except httpx.TimeoutException as exc:
            raise TimeoutError(
                f"POST request timed out for {self.source_name}"
            ) from exc

    def _parse_response(self, response: httpx.Response) -> Any:
        """Parse HTTP response to Python object.

        Args:
            response: httpx Response object.

        Returns:
            Parsed data (JSON dict/list, or text string).
        """
        content_type = response.headers.get("content-type", "")
        if "application/json" in content_type or response.url.path.endswith("/JSON"):
            try:
                return response.json()
            except Exception:
                return response.text
        try:
            return response.json()
        except Exception:
            return response.text

    def _extract_error(self, response: httpx.Response) -> str:
        """Extract error message from an error response.

        Args:
            response: httpx Response object.

        Returns:
            Extracted error message string.
        """
        try:
            data = response.json()
            # PubChem-style error
            fault = data.get("Fault", {})
            message: str = fault.get("Message", "")
            details: list[str] = fault.get("Details", [])
            if details:
                message += " - " + "; ".join(details)
            if message:
                return message
            # Generic message field
            if "message" in data:
                return str(data["message"])
            if "error" in data:
                return str(data["error"])
        except Exception:
            pass
        return response.text[:200] if response.text else f"HTTP {response.status_code}"

    async def close(self) -> None:
        """Close the HTTP client."""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
            self._client = None

    async def __aenter__(self) -> AsyncHTTPClient:
        return self

    async def __aexit__(self, *args: object) -> None:
        await self.close()
