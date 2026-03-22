"""ChemFuse exception hierarchy."""

from __future__ import annotations


class ChemFuseError(Exception):
    """Base exception for all ChemFuse errors."""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(message)


class SourceError(ChemFuseError):
    """Raised when a source adapter encounters an HTTP or API error."""

    def __init__(
        self,
        message: str,
        source: str | None = None,
        url: str | None = None,
        status_code: int | None = None,
    ) -> None:
        self.source = source
        self.url = url
        self.status_code = status_code
        parts = []
        if source:
            parts.append(f"source={source}")
        if url:
            parts.append(f"url={url}")
        if status_code:
            parts.append(f"status={status_code}")
        detail = f" ({', '.join(parts)})" if parts else ""
        super().__init__(f"{message}{detail}")


class NotFoundError(ChemFuseError):
    """Raised when a compound or resource is not found in the source."""

    def __init__(self, message: str, identifier: str | None = None) -> None:
        self.identifier = identifier
        super().__init__(message)


class RateLimitError(ChemFuseError):
    """Raised when a source rate limit is exceeded."""

    def __init__(self, message: str, retry_after: float | None = None) -> None:
        self.retry_after = retry_after
        super().__init__(message)


class ValidationError(ChemFuseError):
    """Raised when input data fails validation."""

    pass


class ChemFuseTimeoutError(ChemFuseError):
    """Raised when a request or operation times out."""

    pass


# Backward-compatible alias; ChemFuseTimeoutError is the canonical name.
TimeoutError = ChemFuseTimeoutError


class OptionalDependencyError(ChemFuseError):
    """Raised when an optional dependency is required but not installed."""

    def __init__(self, package: str, extra: str | None = None) -> None:
        self.package = package
        self.extra = extra
        if extra:
            install_hint = f"pip install chemfuse[{extra}]  # or: pip install {package}"
        else:
            install_hint = f"pip install {package}"
        super().__init__(
            f"Optional dependency '{package}' is not installed. "
            f"Install it with: {install_hint}"
        )
