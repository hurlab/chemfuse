"""Entry point for running the ChemFuse MCP server.

Usage:
    python -m chemfuse.mcp
"""

import asyncio

from chemfuse.mcp.server import main

asyncio.run(main())
