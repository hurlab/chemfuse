"""Prediction models for drug-likeness, molecular property filters, and ADMET."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


class FilterResult(BaseModel):
    """Result of a single drug-likeness filter evaluation."""

    model_config = ConfigDict(populate_by_name=True)

    pass_filter: bool
    violations: list[str]
    details: dict[str, object]


class DrugLikeness(BaseModel):
    """Combined drug-likeness evaluation using all five standard filters.

    Holds FilterResult for each of the five standard drug-likeness filters
    (Lipinski, Veber, Ghose, Egan, Muegge) and an overall_pass flag.
    """

    model_config = ConfigDict(populate_by_name=True)

    lipinski: FilterResult
    veber: FilterResult
    ghose: FilterResult
    egan: FilterResult
    muegge: FilterResult

    @property
    def overall_pass(self) -> bool:
        """True only if all five filters pass."""
        return all(
            [
                self.lipinski.pass_filter,
                self.veber.pass_filter,
                self.ghose.pass_filter,
                self.egan.pass_filter,
                self.muegge.pass_filter,
            ]
        )

    def summary(self) -> dict[str, bool]:
        """Return a dict of filter name to pass/fail boolean.

        Returns:
            Dictionary mapping filter names to pass/fail status.
        """
        return {
            "lipinski": self.lipinski.pass_filter,
            "veber": self.veber.pass_filter,
            "ghose": self.ghose.pass_filter,
            "egan": self.egan.pass_filter,
            "muegge": self.muegge.pass_filter,
        }

    def all_passed(self) -> bool:
        """Return True if all five filters pass.

        Provided as a method for backwards compatibility; use the
        ``overall_pass`` property for idiomatic Pydantic access.
        """
        return self.overall_pass


class PAINSAlert(BaseModel):
    """A single PAINS (Pan Assay Interference Compound) substructure alert."""

    model_config = ConfigDict(populate_by_name=True)

    filter_name: str
    description: str
    matched_atoms: list[int]


class ADMETPrediction(BaseModel):
    """A single ADMET property prediction with confidence and method."""

    model_config = ConfigDict(populate_by_name=True)

    property_name: str
    value: float | str
    unit: str | None = None
    confidence: float | None = None
    # "ml" when admet-ai is used, "rule-based" otherwise
    method: str = "rule-based"
    category: str | None = None  # "high", "medium", "low" or ADMET category


class ADMETProfile(BaseModel):
    """Complete ADMET profile for a single compound."""

    model_config = ConfigDict(populate_by_name=True)

    smiles: str
    predictions: dict[str, ADMETPrediction] = Field(default_factory=dict)
    overall_score: float | None = None
