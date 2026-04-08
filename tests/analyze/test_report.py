"""Tests for analyze/report.py — SAR report generation."""

from __future__ import annotations

import pytest

from chemfuse.models.bioactivity import Bioactivity
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(c1=O)C"


@pytest.fixture
def three_compound_collection() -> CompoundCollection:
    """Collection with aspirin, ibuprofen, and caffeine."""
    aspirin = Compound(
        cid=2244,
        smiles=ASPIRIN_SMILES,
        name="aspirin",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=180.16,
            xlogp=1.2,
            tpsa=63.6,
            hbd_count=1,
            hba_count=4,
            rotatable_bonds=3,
        ),
    )
    ibuprofen = Compound(
        cid=3672,
        smiles=IBUPROFEN_SMILES,
        name="ibuprofen",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=206.28,
            xlogp=3.5,
            tpsa=37.3,
            hbd_count=1,
            hba_count=2,
            rotatable_bonds=4,
        ),
    )
    caffeine = Compound(
        cid=2519,
        smiles=CAFFEINE_SMILES,
        name="caffeine",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=194.19,
            xlogp=-0.07,
            tpsa=58.4,
            hbd_count=0,
            hba_count=6,
            rotatable_bonds=0,
        ),
    )
    return CompoundCollection(
        compounds=[aspirin, ibuprofen, caffeine],
        query="test",
        sources=["pubchem"],
    )


@pytest.fixture
def bioactivity_collection() -> CompoundCollection:
    """Collection with bioactivity data attached to compounds."""
    bio_aspirin = Bioactivity(
        target_name="COX-1",
        activity_type="IC50",
        value=150.0,
        units="nM",
    )
    bio_ibuprofen = Bioactivity(
        target_name="COX-1",
        activity_type="IC50",
        value=800.0,
        units="nM",
    )

    aspirin = Compound(
        cid=2244,
        smiles=ASPIRIN_SMILES,
        name="aspirin",
        sources=["chembl"],
        properties=CompoundProperties(molecular_weight=180.16, xlogp=1.2),
        bioactivities=[bio_aspirin],
    )
    ibuprofen = Compound(
        cid=3672,
        smiles=IBUPROFEN_SMILES,
        name="ibuprofen",
        sources=["chembl"],
        properties=CompoundProperties(molecular_weight=206.28, xlogp=3.5),
        bioactivities=[bio_ibuprofen],
    )
    return CompoundCollection(
        compounds=[aspirin, ibuprofen],
        query="COX-1 inhibitors",
        sources=["chembl"],
    )


# ---------------------------------------------------------------------------
# Tests: generate_sar_report (module-level function)
# ---------------------------------------------------------------------------

class TestGenerateSarReport:
    def test_returns_string(self, three_compound_collection):
        """generate_sar_report returns a non-empty string."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_markdown_format_contains_sections(self, three_compound_collection):
        """Markdown output contains all required section headers."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection, format="markdown")

        assert "## Executive Summary" in result
        assert "## Scaffold Analysis" in result
        assert "## Property Distribution" in result
        assert "## Activity Summary" in result
        assert "## Drug-likeness Profile" in result
        assert "## Key Findings" in result

    def test_markdown_executive_summary_contains_count(self, three_compound_collection):
        """Executive summary reports the correct compound count."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection)
        assert "**Compound count**: 3" in result

    def test_markdown_source_databases_listed(self, three_compound_collection):
        """Executive summary lists source databases."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection)
        assert "pubchem" in result

    def test_text_format_contains_sections(self, three_compound_collection):
        """Plain text output contains all required section headings."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection, format="text")

        assert "EXECUTIVE SUMMARY" in result
        assert "SCAFFOLD ANALYSIS" in result
        assert "PROPERTY DISTRIBUTION" in result
        assert "ACTIVITY SUMMARY" in result
        assert "DRUG-LIKENESS PROFILE" in result
        assert "KEY FINDINGS" in result

    def test_empty_collection_returns_report(self):
        """generate_sar_report handles empty collection gracefully."""
        from chemfuse.analyze.report import generate_sar_report
        empty = CompoundCollection()
        result = generate_sar_report(empty)
        assert isinstance(result, str)
        assert "0" in result or "empty" in result.lower()

    def test_no_activity_data_message(self, three_compound_collection):
        """When no bioactivity data is present, the report notes this."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection)
        # Either the activity summary says "N/A" or the findings mention no data
        assert "N/A" in result or "no bioactivity" in result.lower() or "no ic50" in result.lower()

    def test_with_bioactivity_data_shows_range(self, bioactivity_collection):
        """When bioactivity data is present, the report shows activity range."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(bioactivity_collection, activity_type="IC50")
        # Should show actual values 150 and 800 nM
        assert "150" in result or "800" in result

    def test_with_bioactivity_most_potent_listed(self, bioactivity_collection):
        """Most potent compound (lowest IC50) is listed in the report."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(bioactivity_collection, activity_type="IC50")
        # aspirin is more potent (150 nM < 800 nM)
        assert "aspirin" in result or "150" in result

    def test_target_parameter_in_title(self, bioactivity_collection):
        """When target is specified, it appears in the report header."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(
            bioactivity_collection, activity_type="IC50", target="COX-1"
        )
        assert "COX-1" in result

    def test_collection_method_delegates_to_report(self, three_compound_collection):
        """CompoundCollection.generate_sar_report() delegates to the module function."""
        report = three_compound_collection.generate_sar_report()
        assert isinstance(report, str)
        assert "## Executive Summary" in report

    def test_property_distribution_values(self, three_compound_collection):
        """Property distribution section contains numeric values for MW."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection)
        # Mean MW of aspirin(180.16), ibuprofen(206.28), caffeine(194.19) ≈ 193.54
        assert "180" in result or "193" in result or "206" in result

    def test_key_findings_non_empty(self, three_compound_collection):
        """Key findings section contains at least one bullet point."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection, format="markdown")
        # Find the Key Findings section
        idx = result.find("## Key Findings")
        assert idx >= 0
        section = result[idx:]
        assert "-" in section  # at least one bullet

    def test_text_format_compound_count(self, three_compound_collection):
        """Text format includes compound count in executive summary."""
        from chemfuse.analyze.report import generate_sar_report
        result = generate_sar_report(three_compound_collection, format="text")
        assert "3" in result
