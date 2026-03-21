"""Patent data model for SureChEMBL patent chemistry data."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, field_validator


class Patent(BaseModel):
    """Patent record from SureChEMBL or similar patent chemistry database.

    Represents a single patent document containing chemical compound data.
    """

    model_config = ConfigDict(populate_by_name=True)

    # R-PATENT-02: patent_id must be non-empty
    patent_id: str
    title: str | None = None
    filing_date: str | None = None
    assignee: str | None = None
    jurisdiction: str | None = None
    source: str = "surechembl"

    @field_validator("patent_id")
    @classmethod
    def patent_id_not_empty(cls, v: str) -> str:
        """Validate that patent_id is a non-empty string."""
        if not v or not v.strip():
            raise ValueError("patent_id must be a non-empty string")
        return v

    @property
    def year(self) -> int | None:
        """Extract the filing year from the filing_date string.

        Supports common date formats: YYYY-MM-DD, YYYY/MM/DD, or bare YYYY.

        Returns:
            Integer year, or None if filing_date is absent or cannot be parsed.
        """
        if not self.filing_date:
            return None
        # Extract first 4 digits as year
        stripped = self.filing_date.strip()
        if len(stripped) >= 4 and stripped[:4].isdigit():
            return int(stripped[:4])
        return None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to flat dictionary.

        Returns:
            Dictionary representation of the patent.
        """
        return self.model_dump()


__all__ = ["Patent"]
