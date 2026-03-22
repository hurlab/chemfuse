"""Tests for chemfuse.compute.druglikeness - all five filters and combined evaluation."""

from __future__ import annotations

from chemfuse.models.prediction import DrugLikeness, FilterResult

# ---------------------------------------------------------------------------
# Known compound properties (from PubChem, no RDKit needed)
# ---------------------------------------------------------------------------

ASPIRIN_PROPS = {
    "molecular_weight": 180.16,
    "xlogp": 1.2,
    "tpsa": 63.6,
    "hbd_count": 1,
    "hba_count": 4,
    "rotatable_bonds": 3,
    "heavy_atom_count": 13,
}

CAFFEINE_PROPS = {
    "molecular_weight": 194.19,
    "xlogp": -0.07,
    "tpsa": 58.4,
    "hbd_count": 0,
    "hba_count": 6,
    "rotatable_bonds": 0,
    "heavy_atom_count": 14,
}

# Cyclosporine A: FAIL (MW=1202.61, HBA=23, TPSA=278.8, rotatable=15)
CYCLOSPORINE_PROPS = {
    "molecular_weight": 1202.61,
    "xlogp": 2.9,
    "tpsa": 278.8,
    "hbd_count": 5,
    "hba_count": 23,
    "rotatable_bonds": 15,
    "heavy_atom_count": 88,
}

# A compound that violates MW > 500 but passes everything else (1 Lipinski violation)
ONE_VIOLATION_PROPS = {
    "molecular_weight": 520.0,
    "xlogp": 2.0,
    "tpsa": 80.0,
    "hbd_count": 2,
    "hba_count": 5,
    "rotatable_bonds": 5,
    "heavy_atom_count": 38,
}

# A compound violating two Lipinski rules -> should FAIL Lipinski
TWO_VIOLATIONS_PROPS = {
    "molecular_weight": 520.0,
    "xlogp": 6.0,  # violates LogP <= 5
    "tpsa": 80.0,
    "hbd_count": 2,
    "hba_count": 5,
    "rotatable_bonds": 5,
    "heavy_atom_count": 38,
}


class TestLipinskiFilter:
    """Tests for lipinski_filter()."""

    def test_aspirin_passes(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True
        assert result.violations == []

    def test_caffeine_passes(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(CAFFEINE_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False
        assert len(result.violations) >= 2  # MW and HBA

    def test_one_violation_passes(self):
        """Lipinski allows 1 violation."""
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(ONE_VIOLATION_PROPS)
        assert result.pass_filter is True

    def test_two_violations_fail(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(TWO_VIOLATIONS_PROPS)
        assert result.pass_filter is False
        assert len(result.violations) == 2

    def test_violations_contain_descriptive_messages(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(CYCLOSPORINE_PROPS)
        # Check that MW violation message is informative
        mw_violations = [v for v in result.violations if "MW" in v or "molecular" in v.lower()]
        assert mw_violations, f"Expected MW violation message, got: {result.violations}"
        assert any("1202" in v or "500" in v for v in mw_violations)

    def test_returns_filter_result_instance(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(ASPIRIN_PROPS)
        assert isinstance(result, FilterResult)

    def test_details_contains_property_info(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter(ASPIRIN_PROPS)
        assert "molecular_weight" in result.details
        assert "logp" in result.details

    def test_missing_properties_returns_fail(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        result = lipinski_filter({})
        assert result.pass_filter is False
        assert len(result.violations) > 0

    def test_boundary_mw_exactly_500_passes(self):
        from chemfuse.compute.druglikeness import lipinski_filter
        props = {**ASPIRIN_PROPS, "molecular_weight": 500.0}
        result = lipinski_filter(props)
        # MW == 500 should pass (threshold is <=500)
        mw_violations = [v for v in result.violations if "MW" in v or "500" in v]
        assert not mw_violations


class TestVeberFilter:
    """Tests for veber_filter()."""

    def test_aspirin_passes(self):
        from chemfuse.compute.druglikeness import veber_filter
        result = veber_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_caffeine_passes(self):
        from chemfuse.compute.druglikeness import veber_filter
        result = veber_filter(CAFFEINE_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        from chemfuse.compute.druglikeness import veber_filter
        result = veber_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False

    def test_high_tpsa_fails(self):
        from chemfuse.compute.druglikeness import veber_filter
        props = {"tpsa": 150.0, "rotatable_bonds": 5}
        result = veber_filter(props)
        assert result.pass_filter is False
        assert any("TPSA" in v for v in result.violations)

    def test_high_rotatable_bonds_fails(self):
        from chemfuse.compute.druglikeness import veber_filter
        props = {"tpsa": 80.0, "rotatable_bonds": 11}
        result = veber_filter(props)
        assert result.pass_filter is False
        assert any("RotBonds" in v or "rotatable" in v.lower() for v in result.violations)

    def test_boundary_tpsa_140_passes(self):
        from chemfuse.compute.druglikeness import veber_filter
        props = {"tpsa": 140.0, "rotatable_bonds": 5}
        result = veber_filter(props)
        assert result.pass_filter is True

    def test_boundary_rotatable_bonds_10_passes(self):
        from chemfuse.compute.druglikeness import veber_filter
        props = {"tpsa": 80.0, "rotatable_bonds": 10}
        result = veber_filter(props)
        assert result.pass_filter is True


class TestGhoseFilter:
    """Tests for ghose_filter()."""

    def test_aspirin_passes(self):
        from chemfuse.compute.druglikeness import ghose_filter
        result = ghose_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_caffeine_passes(self):
        from chemfuse.compute.druglikeness import ghose_filter
        result = ghose_filter(CAFFEINE_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails_mw(self):
        from chemfuse.compute.druglikeness import ghose_filter
        result = ghose_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False

    def test_too_small_mw_fails(self):
        from chemfuse.compute.druglikeness import ghose_filter
        props = {"molecular_weight": 100.0, "xlogp": 1.0, "heavy_atom_count": 30}
        result = ghose_filter(props)
        assert result.pass_filter is False

    def test_too_large_mw_fails(self):
        from chemfuse.compute.druglikeness import ghose_filter
        props = {"molecular_weight": 500.0, "xlogp": 1.0, "heavy_atom_count": 30}
        result = ghose_filter(props)
        assert result.pass_filter is False

    def test_out_of_logp_range_fails(self):
        from chemfuse.compute.druglikeness import ghose_filter
        props = {**ASPIRIN_PROPS, "xlogp": 6.0}
        result = ghose_filter(props)
        assert result.pass_filter is False

    def test_boundary_mw_160_passes(self):
        from chemfuse.compute.druglikeness import ghose_filter
        props = {
            "molecular_weight": 160.0,
            "xlogp": 1.0,
            "heavy_atom_count": 30,
        }
        result = ghose_filter(props)
        # Other violations may apply (atom count), but MW itself should pass
        mw_violations = [v for v in result.violations if "MW" in v]
        assert not mw_violations


class TestEganFilter:
    """Tests for egan_filter()."""

    def test_aspirin_passes(self):
        from chemfuse.compute.druglikeness import egan_filter
        result = egan_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_caffeine_passes(self):
        from chemfuse.compute.druglikeness import egan_filter
        result = egan_filter(CAFFEINE_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails_tpsa(self):
        from chemfuse.compute.druglikeness import egan_filter
        result = egan_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False
        assert any("TPSA" in v for v in result.violations)

    def test_high_logp_fails(self):
        from chemfuse.compute.druglikeness import egan_filter
        props = {"xlogp": 6.0, "tpsa": 80.0}
        result = egan_filter(props)
        assert result.pass_filter is False
        assert any("LogP" in v for v in result.violations)

    def test_boundary_tpsa_131_6_passes(self):
        from chemfuse.compute.druglikeness import egan_filter
        props = {"xlogp": 1.0, "tpsa": 131.6}
        result = egan_filter(props)
        assert result.pass_filter is True

    def test_negative_logp_in_range_passes(self):
        from chemfuse.compute.druglikeness import egan_filter
        props = {"xlogp": -0.5, "tpsa": 60.0}
        result = egan_filter(props)
        assert result.pass_filter is True

    def test_logp_too_negative_fails(self):
        from chemfuse.compute.druglikeness import egan_filter
        props = {"xlogp": -2.0, "tpsa": 60.0}
        result = egan_filter(props)
        assert result.pass_filter is False


class TestMueggeFilter:
    """Tests for muegge_filter()."""

    def test_aspirin_passes(self):
        from chemfuse.compute.druglikeness import muegge_filter
        # Aspirin needs ring and heteroatom data. Use SMILES for RDKit computation.
        result = muegge_filter(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert result.pass_filter is True

    def test_caffeine_passes(self):
        from chemfuse.compute.druglikeness import muegge_filter
        result = muegge_filter(CAFFEINE_PROPS, smiles="Cn1c(=O)c2c(ncn2C)n(c1=O)C")
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        from chemfuse.compute.druglikeness import muegge_filter
        result = muegge_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False

    def test_mw_too_small_fails(self):
        from chemfuse.compute.druglikeness import muegge_filter
        props = {**ASPIRIN_PROPS, "molecular_weight": 150.0}
        result = muegge_filter(props)
        assert result.pass_filter is False

    def test_mw_too_large_fails(self):
        from chemfuse.compute.druglikeness import muegge_filter
        props = {**ASPIRIN_PROPS, "molecular_weight": 650.0}
        result = muegge_filter(props)
        assert result.pass_filter is False

    def test_returns_filter_result(self):
        from chemfuse.compute.druglikeness import muegge_filter
        result = muegge_filter(ASPIRIN_PROPS)
        assert isinstance(result, FilterResult)


class TestCheckDrugLikeness:
    """Tests for check_drug_likeness() combined evaluation."""

    def test_aspirin_returns_drug_likeness(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(result, DrugLikeness)

    def test_all_five_filters_present(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(result.lipinski, FilterResult)
        assert isinstance(result.veber, FilterResult)
        assert isinstance(result.ghose, FilterResult)
        assert isinstance(result.egan, FilterResult)
        assert isinstance(result.muegge, FilterResult)

    def test_aspirin_all_pass(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert result.lipinski.pass_filter is True
        assert result.veber.pass_filter is True

    def test_overall_pass_property(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        # overall_pass is True only if all filters pass
        expected = all([
            result.lipinski.pass_filter,
            result.veber.pass_filter,
            result.ghose.pass_filter,
            result.egan.pass_filter,
            result.muegge.pass_filter,
        ])
        assert result.overall_pass == expected

    def test_summary_returns_dict_of_bools(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        summary = result.summary()
        assert isinstance(summary, dict)
        # Five core filters are always present; pains is present when SMILES given.
        core_keys = {"lipinski", "veber", "ghose", "egan", "muegge"}
        assert core_keys.issubset(summary.keys())
        for key, val in summary.items():
            assert isinstance(val, bool), f"summary[{key!r}] is not bool"

    def test_cyclosporine_fails_overall(self):
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(CYCLOSPORINE_PROPS)
        assert result.overall_pass is False

    def test_without_rdkit_uses_provided_props(self):
        """Filters must work using only provided properties (no RDKit required)."""
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import check_drug_likeness
        original = dl_mod._RDKIT_AVAILABLE
        dl_mod._RDKIT_AVAILABLE = False
        try:
            result = check_drug_likeness(ASPIRIN_PROPS)
            assert isinstance(result, DrugLikeness)
            # Aspirin passes Lipinski and Veber even without RDKit
            assert result.lipinski.pass_filter is True
            assert result.veber.pass_filter is True
        finally:
            dl_mod._RDKIT_AVAILABLE = original

    def test_with_smiles_populates_pains_and_qed(self):
        """check_drug_likeness populates pains and qed when smiles is given."""
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert result.pains is not None, "pains should be populated when SMILES provided"
        assert result.qed is not None, "qed should be populated when SMILES provided"

    def test_without_smiles_pains_and_qed_are_none(self):
        """check_drug_likeness omits pains and qed when no SMILES is provided."""
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS)
        assert result.pains is None
        assert result.qed is None

    def test_summary_includes_pains_when_available(self):
        """summary() includes 'pains' key when PAINS result is present."""
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        summary = result.summary()
        assert "pains" in summary
        assert isinstance(summary["pains"], bool)

    def test_summary_excludes_pains_when_unavailable(self):
        """summary() omits 'pains' key when no SMILES was provided."""
        from chemfuse.compute.druglikeness import check_drug_likeness
        result = check_drug_likeness(ASPIRIN_PROPS)
        summary = result.summary()
        assert "pains" not in summary


class TestPAINSFilter:
    """Tests for pains_filter()."""

    def test_aspirin_passes(self):
        """Aspirin has no PAINS alerts."""
        from chemfuse.compute.druglikeness import pains_filter
        result = pains_filter("CC(=O)Oc1ccccc1C(=O)O")
        assert result.pass_filter is True
        assert result.violations == []
        assert result.details["pains_match"] is None

    def test_rhodanine_fails(self):
        """Rhodanine is a well-known PAINS scaffold."""
        from chemfuse.compute.druglikeness import pains_filter
        result = pains_filter("O=C1CSC(=S)N1")
        assert result.pass_filter is False
        assert len(result.violations) == 1
        assert result.details["pains_match"] is not None

    def test_catechol_fails(self):
        """Catechol is a well-known PAINS scaffold (catechol_A alert)."""
        from chemfuse.compute.druglikeness import pains_filter
        result = pains_filter("Oc1ccccc1O")
        assert result.pass_filter is False

    def test_returns_filter_result(self):
        """pains_filter returns a FilterResult instance."""
        from chemfuse.compute.druglikeness import pains_filter
        from chemfuse.models.prediction import FilterResult
        result = pains_filter("CC(=O)Oc1ccccc1C(=O)O")
        assert isinstance(result, FilterResult)

    def test_violation_message_contains_pattern_name(self):
        """The violation string names the matched PAINS pattern."""
        from chemfuse.compute.druglikeness import pains_filter
        result = pains_filter("O=C1CSC(=S)N1")
        assert result.violations, "Expected at least one violation message"
        assert result.details["pains_match"] in result.violations[0]

    def test_invalid_smiles_raises_value_error(self):
        """pains_filter raises ValueError for unparseable SMILES."""
        import pytest

        from chemfuse.compute.druglikeness import pains_filter
        with pytest.raises(ValueError, match="Could not parse SMILES"):
            pains_filter("not_a_smiles!!!")

    def test_no_rdkit_raises_import_error(self):
        """pains_filter raises ImportError when RDKit is not available."""
        import pytest

        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import pains_filter
        original = dl_mod._RDKIT_AVAILABLE
        dl_mod._RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="RDKit is required"):
                pains_filter("CC(=O)Oc1ccccc1C(=O)O")
        finally:
            dl_mod._RDKIT_AVAILABLE = original


class TestQEDScore:
    """Tests for qed_score()."""

    def test_aspirin_moderate_qed(self):
        """Aspirin has a moderate QED (>0.4)."""
        from chemfuse.compute.druglikeness import qed_score
        result = qed_score("CC(=O)Oc1ccccc1C(=O)O")
        assert result["qed"] > 0.4, f"Expected QED > 0.4, got {result['qed']}"

    def test_benzene_low_qed(self):
        """Benzene is not drug-like; QED should be classified low or medium."""
        from chemfuse.compute.druglikeness import qed_score
        result = qed_score("c1ccccc1")
        # Benzene QED is empirically low (~0.17)
        assert result["qed"] < 0.5, f"Expected benzene QED < 0.5, got {result['qed']}"

    def test_ibuprofen_high_qed(self):
        """Ibuprofen is a well-known oral drug; QED should be high."""
        from chemfuse.compute.druglikeness import qed_score
        # Ibuprofen: MW 206, LogP 3.97, well-studied QED ~0.73
        result = qed_score("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        assert result["qed"] > 0.6, f"Expected ibuprofen QED > 0.6, got {result['qed']}"

    def test_classification_high(self):
        """QED >= 0.67 is classified 'high'."""
        from chemfuse.compute.druglikeness import qed_score
        result = qed_score("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        if result["qed"] >= 0.67:
            assert result["classification"] == "high"

    def test_classification_low(self):
        """QED < 0.33 is classified 'low'."""
        from chemfuse.compute.druglikeness import qed_score
        result = qed_score("c1ccccc1")
        if result["qed"] < 0.33:
            assert result["classification"] == "low"

    def test_returns_dict_with_expected_keys(self):
        """qed_score returns a dict with 'qed' and 'classification' keys."""
        from chemfuse.compute.druglikeness import qed_score
        result = qed_score("CC(=O)Oc1ccccc1C(=O)O")
        assert "qed" in result
        assert "classification" in result
        assert isinstance(result["qed"], float)
        assert result["classification"] in ("high", "medium", "low")

    def test_qed_in_unit_interval(self):
        """QED score is always between 0 and 1 inclusive."""
        from chemfuse.compute.druglikeness import qed_score
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
            "c1ccccc1",               # benzene
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
        ]:
            result = qed_score(smi)
            assert 0.0 <= result["qed"] <= 1.0, f"QED out of range for {smi}: {result['qed']}"

    def test_invalid_smiles_raises_value_error(self):
        """qed_score raises ValueError for unparseable SMILES."""
        import pytest

        from chemfuse.compute.druglikeness import qed_score
        with pytest.raises(ValueError, match="Could not parse SMILES"):
            qed_score("not_a_smiles!!!")

    def test_no_rdkit_raises_import_error(self):
        """qed_score raises ImportError when RDKit is not available."""
        import pytest

        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import qed_score
        original = dl_mod._RDKIT_AVAILABLE
        dl_mod._RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="RDKit is required"):
                qed_score("CC(=O)Oc1ccccc1C(=O)O")
        finally:
            dl_mod._RDKIT_AVAILABLE = original
