"""TargetAssociation data model for Open Targets disease-target associations."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


class TargetAssociation(BaseModel):
    """Disease-target association from Open Targets Platform.

    Represents a single association between a drug target (gene/protein)
    and a disease, including the association score and evidence.
    """

    model_config = ConfigDict(populate_by_name=True)

    target_id: str
    target_name: str
    disease_id: str | None = None
    disease_name: str | None = None
    # R-TARGET-02: association_score must be between 0.0 and 1.0
    association_score: float | None = Field(default=None, ge=0.0, le=1.0)
    evidence_count: int | None = None
    evidence_types: list[str] = Field(default_factory=list)
    source: str = "opentargets"

    @field_validator("association_score", mode="before")
    @classmethod
    def clamp_score(cls, v: Any) -> Any:
        """Clamp association_score to [0.0, 1.0] range before validation."""
        if v is None:
            return v
        try:
            fv = float(v)
            return max(0.0, min(1.0, fv))
        except (TypeError, ValueError):
            return v

    @property
    def strength(self) -> str:
        """Categorize association strength based on association_score.

        Returns:
            'strong' if score >= 0.7,
            'moderate' if score >= 0.4,
            'weak' if score < 0.4,
            or 'unknown' if score is None.
        """
        if self.association_score is None:
            return "unknown"
        if self.association_score >= 0.7:
            return "strong"
        if self.association_score >= 0.4:
            return "moderate"
        return "weak"

    def to_dict(self) -> dict[str, Any]:
        """Serialize to flat dictionary.

        Returns:
            Dictionary representation of the target association.
        """
        return self.model_dump()


__all__ = ["TargetAssociation"]
