"""SourceAdapter abstract base class and protocol."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from chemfuse.models.compound import Compound


class SourceAdapter(ABC):
    """Abstract base class for chemical database source adapters.

    All source adapters must implement this interface to be compatible
    with the ChemFuse source registry and search framework.
    """

    # Must be set by subclasses
    name: str = ""
    base_url: str = ""
    rate_limit: float = 1.0  # requests per second

    @abstractmethod
    async def search(
        self,
        query: str,
        query_type: str = "name",
    ) -> list[Compound]:
        """Search for compounds by query.

        Args:
            query: Search query string.
            query_type: Type of query - one of 'name', 'smiles', 'formula',
                        'identifier', 'inchi', 'cid'.

        Returns:
            List of Compound objects matching the query.
        """
        ...

    @abstractmethod
    async def get_by_id(self, identifier: str) -> Compound | None:
        """Get a compound by its source-specific identifier.

        Args:
            identifier: Source-specific identifier (e.g., CID for PubChem).

        Returns:
            Compound object, or None if not found.
        """
        ...

    @abstractmethod
    async def get_properties(self, identifier: str) -> dict:
        """Get detailed properties for a compound.

        Args:
            identifier: Source-specific identifier.

        Returns:
            Dictionary of property names to values.
        """
        ...

    @abstractmethod
    def is_available(self) -> bool:
        """Check if the source API is currently available.

        Returns:
            True if the source is reachable and functional.
        """
        ...
