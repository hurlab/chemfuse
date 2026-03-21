"""Tests for TargetAssociation model."""

from __future__ import annotations

import pytest

from chemfuse.models.target import TargetAssociation


class TestTargetAssociationCreation:
    """Tests for TargetAssociation model creation."""

    def test_create_minimal(self) -> None:
        """Create TargetAssociation with required fields only."""
        assoc = TargetAssociation(target_id="ENSG00000095303", target_name="PTGS1")
        assert assoc.target_id == "ENSG00000095303"
        assert assoc.target_name == "PTGS1"
        assert assoc.disease_id is None
        assert assoc.disease_name is None
        assert assoc.association_score is None
        assert assoc.evidence_count is None
        assert assoc.evidence_types == []
        assert assoc.source == "opentargets"

    def test_create_full(self) -> None:
        """Create TargetAssociation with all fields."""
        assoc = TargetAssociation(
            target_id="ENSG00000095303",
            target_name="PTGS1",
            disease_id="EFO_0003785",
            disease_name="pain",
            association_score=0.85,
            evidence_count=42,
            evidence_types=["genetic_association", "somatic_mutation"],
            source="opentargets",
        )
        assert assoc.target_id == "ENSG00000095303"
        assert assoc.target_name == "PTGS1"
        assert assoc.disease_id == "EFO_0003785"
        assert assoc.disease_name == "pain"
        assert assoc.association_score == pytest.approx(0.85)
        assert assoc.evidence_count == 42
        assert assoc.evidence_types == ["genetic_association", "somatic_mutation"]

    def test_default_source(self) -> None:
        """Default source is 'opentargets'."""
        assoc = TargetAssociation(target_id="T1", target_name="Gene1")
        assert assoc.source == "opentargets"

    def test_custom_source(self) -> None:
        """Custom source is preserved."""
        assoc = TargetAssociation(target_id="T1", target_name="Gene1", source="custom_db")
        assert assoc.source == "custom_db"


class TestTargetAssociationValidation:
    """Tests for TargetAssociation field validation."""

    def test_score_zero_is_valid(self) -> None:
        """association_score=0.0 is valid."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.0
        )
        assert assoc.association_score == pytest.approx(0.0)

    def test_score_one_is_valid(self) -> None:
        """association_score=1.0 is valid."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=1.0
        )
        assert assoc.association_score == pytest.approx(1.0)

    def test_score_clamped_above_one(self) -> None:
        """Score above 1.0 is clamped to 1.0."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=1.5
        )
        assert assoc.association_score == pytest.approx(1.0)

    def test_score_clamped_below_zero(self) -> None:
        """Score below 0.0 is clamped to 0.0."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=-0.1
        )
        assert assoc.association_score == pytest.approx(0.0)

    def test_score_none_is_valid(self) -> None:
        """None score is valid."""
        assoc = TargetAssociation(target_id="T1", target_name="G1", association_score=None)
        assert assoc.association_score is None

    def test_score_mid_range(self) -> None:
        """Mid-range score is stored unchanged."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.5
        )
        assert assoc.association_score == pytest.approx(0.5)


class TestTargetAssociationStrengthProperty:
    """Tests for TargetAssociation.strength property."""

    def test_strong_score_0_8(self) -> None:
        """Score 0.8 is classified as 'strong' (>= 0.7)."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.8
        )
        assert assoc.strength == "strong"

    def test_strong_exact_boundary(self) -> None:
        """Score exactly 0.7 is classified as 'strong'."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.7
        )
        assert assoc.strength == "strong"

    def test_moderate_score_0_5(self) -> None:
        """Score 0.5 is classified as 'moderate' (>= 0.4)."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.5
        )
        assert assoc.strength == "moderate"

    def test_moderate_exact_boundary(self) -> None:
        """Score exactly 0.4 is classified as 'moderate'."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.4
        )
        assert assoc.strength == "moderate"

    def test_weak_score_0_2(self) -> None:
        """Score 0.2 is classified as 'weak' (< 0.4)."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.2
        )
        assert assoc.strength == "weak"

    def test_weak_score_0_0(self) -> None:
        """Score 0.0 is classified as 'weak'."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.0
        )
        assert assoc.strength == "weak"

    def test_none_score_returns_unknown(self) -> None:
        """None score returns 'unknown'."""
        assoc = TargetAssociation(target_id="T1", target_name="G1", association_score=None)
        assert assoc.strength == "unknown"

    def test_score_just_below_strong(self) -> None:
        """Score 0.699 is 'moderate', not 'strong'."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.699
        )
        assert assoc.strength == "moderate"

    def test_score_just_below_moderate(self) -> None:
        """Score 0.399 is 'weak', not 'moderate'."""
        assoc = TargetAssociation(
            target_id="T1", target_name="G1", association_score=0.399
        )
        assert assoc.strength == "weak"


class TestTargetAssociationSerialization:
    """Tests for TargetAssociation serialization."""

    def test_to_dict_returns_dict(self) -> None:
        """to_dict returns a dictionary."""
        assoc = TargetAssociation(
            target_id="T1",
            target_name="Gene1",
            disease_name="condition",
            association_score=0.75,
        )
        d = assoc.to_dict()
        assert isinstance(d, dict)

    def test_to_dict_contains_all_fields(self) -> None:
        """to_dict contains all model fields."""
        assoc = TargetAssociation(
            target_id="T1",
            target_name="Gene1",
            disease_id="EFO_001",
            disease_name="pain",
            association_score=0.75,
            evidence_count=10,
        )
        d = assoc.to_dict()
        assert d["target_id"] == "T1"
        assert d["target_name"] == "Gene1"
        assert d["disease_id"] == "EFO_001"
        assert d["disease_name"] == "pain"
        assert d["association_score"] == pytest.approx(0.75)
        assert d["evidence_count"] == 10

    def test_model_dump_matches_to_dict(self) -> None:
        """model_dump output matches to_dict."""
        assoc = TargetAssociation(target_id="T1", target_name="Gene1", association_score=0.5)
        assert assoc.to_dict() == assoc.model_dump()

    def test_evidence_types_list_in_dict(self) -> None:
        """evidence_types list is preserved in serialization."""
        assoc = TargetAssociation(
            target_id="T1",
            target_name="Gene1",
            evidence_types=["genetic_association"],
        )
        d = assoc.to_dict()
        assert d["evidence_types"] == ["genetic_association"]
