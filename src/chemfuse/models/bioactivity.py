"""Bioactivity and BindingMeasurement data models for ChemFuse."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, field_validator

# Recognized activity types per SPEC R-BIOACT-02
RECOGNIZED_ACTIVITY_TYPES = frozenset({
    "IC50", "Ki", "Kd", "EC50", "ED50", "GI50", "other",
})


class Bioactivity(BaseModel):
    """Bioactivity measurement from ChEMBL or other sources.

    Represents a single bioactivity data point linking a compound to a
    biological target with a measured activity value.
    """

    model_config = ConfigDict(populate_by_name=True)

    # Target information
    target_name: str
    target_id: str | None = None  # e.g., ChEMBL target ID

    # Activity measurement
    activity_type: str  # IC50, Ki, Kd, EC50, ED50, GI50, or other
    value: float | None = None
    units: str | None = None
    relation: str | None = None  # =, <, >, <=, >=
    assay_type: str | None = None  # Binding, Functional, ADME

    # Source metadata
    source: str = ""  # e.g., "chembl"
    reference: str | None = None  # PubMed ID or URL

    @field_validator("activity_type", mode="before")
    @classmethod
    def normalize_activity_type(cls, v: str) -> str:
        """Normalize activity type to recognized values or 'other'."""
        if v.upper() in {t.upper() for t in RECOGNIZED_ACTIVITY_TYPES}:
            # Return the canonical capitalization
            for canonical in RECOGNIZED_ACTIVITY_TYPES:
                if canonical.upper() == v.upper():
                    return canonical
        return "other"

    def to_dict(self) -> dict:
        """Serialize to flat dictionary.

        Returns:
            Dictionary representation of the bioactivity.
        """
        return self.model_dump()


class BindingMeasurement(BaseModel):
    """Binding affinity measurement from BindingDB or similar sources.

    Aggregates Ki, Kd, IC50, and EC50 values for a single compound-target
    interaction into a unified record.
    """

    model_config = ConfigDict(populate_by_name=True)

    # Target information
    target_name: str
    target_uniprot: str | None = None  # UniProt accession

    # Binding affinity values (all in nM unless specified by units)
    ki: float | None = None
    kd: float | None = None
    ic50: float | None = None
    ec50: float | None = None

    # Units (if non-standard)
    ki_units: str | None = None
    kd_units: str | None = None
    ic50_units: str | None = None
    ec50_units: str | None = None

    # Inequality flags (e.g., ">10000" stores value=10000, inequality=">")
    ki_relation: str | None = None
    kd_relation: str | None = None
    ic50_relation: str | None = None
    ec50_relation: str | None = None

    # Source metadata
    source: str = "bindingdb"
    reference: str | None = None

    @property
    def best_affinity(self) -> float | None:
        """Return the strongest (lowest) available binding affinity value.

        Compares Ki, Kd, IC50, and EC50 and returns the minimum non-None value.
        Lower values indicate stronger binding.

        Returns:
            Minimum affinity value, or None if no measurements are available.
        """
        values = [v for v in [self.ki, self.kd, self.ic50, self.ec50] if v is not None]
        return min(values) if values else None

    def to_dict(self) -> dict:
        """Serialize to flat dictionary.

        Returns:
            Dictionary representation of the binding measurement.
        """
        return self.model_dump()


__all__ = [
    "Bioactivity",
    "BindingMeasurement",
    "RECOGNIZED_ACTIVITY_TYPES",
]
