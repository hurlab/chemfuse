# Builder stage: install dependencies in a separate layer for caching
FROM python:3.12-slim AS builder

WORKDIR /build

# Install build dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        libxrender1 \
        libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Copy package metadata for dependency installation
COPY pyproject.toml README.md LICENSE ./

# Create a minimal src layout so pip can resolve the package
COPY src/ src/

# Install the full package with all optional dependencies
RUN pip install --no-cache-dir --prefix=/install ".[all]"

# Runtime stage: minimal production image
FROM python:3.12-slim

WORKDIR /app

# Install runtime system libraries required by RDKit
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN groupadd -r appuser && useradd -r -g appuser -d /app -s /sbin/nologin appuser

# Copy installed packages from builder stage
COPY --from=builder /install /usr/local

# Copy application source
COPY src/ src/

# Set ownership to non-root user
RUN chown -R appuser:appuser /app

# Create cache directory
RUN mkdir -p /app/.chemfuse && chown appuser:appuser /app/.chemfuse

USER appuser

# Expose Streamlit default port
EXPOSE 8501

# Health check: verify chemfuse can be imported and the web server responds
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Default command: start the web UI
CMD ["chemfuse", "web", "--host", "0.0.0.0"]
