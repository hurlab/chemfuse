"""Tests for Patent model."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemfuse.models.patent import Patent


class TestPatentCreation:
    """Tests for Patent model creation."""

    def test_create_minimal(self) -> None:
        """Create Patent with required fields only."""
        patent = Patent(patent_id="US-4568679-A")
        assert patent.patent_id == "US-4568679-A"
        assert patent.title is None
        assert patent.filing_date is None
        assert patent.assignee is None
        assert patent.jurisdiction is None
        assert patent.source == "surechembl"

    def test_create_full(self) -> None:
        """Create Patent with all fields."""
        patent = Patent(
            patent_id="US-4568679-A",
            title="Coated aspirin tablet",
            filing_date="1984-03-15",
            assignee="Bayer AG",
            jurisdiction="US",
            source="surechembl",
        )
        assert patent.patent_id == "US-4568679-A"
        assert patent.title == "Coated aspirin tablet"
        assert patent.filing_date == "1984-03-15"
        assert patent.assignee == "Bayer AG"
        assert patent.jurisdiction == "US"

    def test_default_source(self) -> None:
        """Default source is 'surechembl'."""
        patent = Patent(patent_id="EP123")
        assert patent.source == "surechembl"

    def test_custom_source(self) -> None:
        """Custom source is preserved."""
        patent = Patent(patent_id="EP123", source="custom_db")
        assert patent.source == "custom_db"


class TestPatentValidation:
    """Tests for Patent field validation (R-PATENT-02)."""

    def test_empty_patent_id_raises(self) -> None:
        """Empty patent_id raises ValidationError."""
        with pytest.raises(ValidationError):
            Patent(patent_id="")

    def test_whitespace_only_patent_id_raises(self) -> None:
        """Whitespace-only patent_id raises ValidationError."""
        with pytest.raises(ValidationError):
            Patent(patent_id="   ")

    def test_valid_patent_id_with_spaces(self) -> None:
        """Patent ID with non-whitespace content is valid."""
        patent = Patent(patent_id="US 4568679 A")
        assert patent.patent_id == "US 4568679 A"

    def test_numeric_patent_id_as_string(self) -> None:
        """Numeric-looking patent ID strings are valid."""
        patent = Patent(patent_id="123456")
        assert patent.patent_id == "123456"


class TestPatentYearProperty:
    """Tests for Patent.year property (R-PATENT-04)."""

    def test_iso_date_extracts_year(self) -> None:
        """YYYY-MM-DD format returns correct year."""
        patent = Patent(patent_id="US-1", filing_date="1984-03-15")
        assert patent.year == 1984

    def test_slash_date_extracts_year(self) -> None:
        """YYYY/MM/DD format returns correct year."""
        patent = Patent(patent_id="US-1", filing_date="2001/05/20")
        assert patent.year == 2001

    def test_year_only_string(self) -> None:
        """Bare YYYY string returns correct year."""
        patent = Patent(patent_id="US-1", filing_date="1999")
        assert patent.year == 1999

    def test_none_filing_date_returns_none(self) -> None:
        """None filing_date returns None year."""
        patent = Patent(patent_id="US-1")
        assert patent.year is None

    def test_empty_filing_date_returns_none(self) -> None:
        """Empty filing_date cannot be parsed; returns None if not parseable."""
        patent = Patent(patent_id="US-1", filing_date="")
        assert patent.year is None

    def test_non_numeric_start_returns_none(self) -> None:
        """Filing date that doesn't start with 4 digits returns None."""
        patent = Patent(patent_id="US-1", filing_date="March 1984")
        assert patent.year is None

    def test_leading_spaces_stripped(self) -> None:
        """Leading spaces in filing_date are stripped before parsing."""
        patent = Patent(patent_id="US-1", filing_date="  2010-06-01")
        assert patent.year == 2010

    def test_year_2024(self) -> None:
        """Recent year is parsed correctly."""
        patent = Patent(patent_id="US-1", filing_date="2024-12-31")
        assert patent.year == 2024


class TestPatentSerialization:
    """Tests for Patent serialization (R-PATENT-03)."""

    def test_to_dict_returns_dict(self) -> None:
        """to_dict returns a dictionary."""
        patent = Patent(
            patent_id="US-4568679-A",
            title="Aspirin coating",
            filing_date="1984-03-15",
        )
        d = patent.to_dict()
        assert isinstance(d, dict)

    def test_to_dict_contains_all_fields(self) -> None:
        """to_dict contains all model fields."""
        patent = Patent(
            patent_id="US-4568679-A",
            title="Coated aspirin tablet",
            filing_date="1984-03-15",
            assignee="Bayer AG",
            jurisdiction="US",
        )
        d = patent.to_dict()
        assert d["patent_id"] == "US-4568679-A"
        assert d["title"] == "Coated aspirin tablet"
        assert d["filing_date"] == "1984-03-15"
        assert d["assignee"] == "Bayer AG"
        assert d["jurisdiction"] == "US"

    def test_model_dump_matches_to_dict(self) -> None:
        """model_dump output matches to_dict."""
        patent = Patent(patent_id="US-123", filing_date="2000-01-01")
        assert patent.to_dict() == patent.model_dump()

    def test_optional_fields_none_in_dict(self) -> None:
        """Optional fields are None in dict when not provided."""
        patent = Patent(patent_id="US-1")
        d = patent.to_dict()
        assert d["title"] is None
        assert d["filing_date"] is None
        assert d["assignee"] is None
        assert d["jurisdiction"] is None
