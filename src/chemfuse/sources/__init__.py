"""ChemFuse source adapter registry."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from chemfuse.sources._base import SourceAdapter


class SourceRegistry:
    """Registry for source adapters with lazy loading.

    Maps source names to adapter classes or instances. Adapters are
    instantiated lazily on first access.
    """

    def __init__(self) -> None:
        self._adapters: dict[str, SourceAdapter] = {}
        self._adapter_classes: dict[str, type] = {}

    def register(self, name: str, adapter_class: type) -> None:
        """Register a source adapter class.

        Args:
            name: Source name identifier.
            adapter_class: SourceAdapter subclass to register.
        """
        self._adapter_classes[name] = adapter_class
        # Clear cached instance if re-registering
        if name in self._adapters:
            del self._adapters[name]

    def get(self, name: str) -> SourceAdapter:
        """Get a source adapter instance by name.

        Instantiates the adapter on first access (lazy loading).

        Args:
            name: Source name identifier.

        Returns:
            Instantiated SourceAdapter.

        Raises:
            KeyError: If the source name is not registered.
        """
        if name not in self._adapters:
            if name not in self._adapter_classes:
                available = ", ".join(sorted(self._adapter_classes.keys()))
                raise KeyError(
                    f"Source '{name}' not registered. Available sources: {available}"
                )
            self._adapters[name] = self._adapter_classes[name]()
        return self._adapters[name]

    def list(self) -> list[str]:
        """List all registered source names.

        Returns:
            Sorted list of registered source names.
        """
        return sorted(self._adapter_classes.keys())

    def is_registered(self, name: str) -> bool:
        """Check if a source name is registered.

        Args:
            name: Source name to check.

        Returns:
            True if the source is registered.
        """
        return name in self._adapter_classes


# Global registry singleton
registry = SourceRegistry()


def _register_defaults() -> None:
    """Register default source adapters."""
    from chemfuse.sources.bindingdb import BindingDBAdapter
    from chemfuse.sources.chembl import ChEMBLAdapter
    from chemfuse.sources.opentargets import OpenTargetsAdapter
    from chemfuse.sources.pubchem import PubChemAdapter
    from chemfuse.sources.surechembl import SureChEMBLAdapter
    from chemfuse.sources.unichem import UniChemAdapter

    registry.register("pubchem", PubChemAdapter)
    registry.register("chembl", ChEMBLAdapter)
    registry.register("unichem", UniChemAdapter)
    registry.register("bindingdb", BindingDBAdapter)
    registry.register("opentargets", OpenTargetsAdapter)
    registry.register("surechembl", SureChEMBLAdapter)


_register_defaults()


__all__ = ["SourceRegistry", "registry"]
