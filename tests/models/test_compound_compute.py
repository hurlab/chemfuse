"""Tests for compute methods on Compound and CompoundCollection (SPEC-CF-002)."""

from __future__ import annotations

import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(c1=O)C"
CYCLOSPORINE_SMILES = (
    "CC[C@@H]1OC(=O)[C@H](CC(=O)c2ccccc2)NC(=O)[C@@H](CC(C)C)N(C)C(=O)[C@H]"
    "(Cc2ccccc2)NC(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](Cc2cccnc2)NC(=O)[C@@H]"
    "(CCCC)NC(=O)CN(C)C1=O"
)


def _make_aspirin() -> Compound:
    return Compound(
        cid=2244,
        smiles=ASPIRIN_SMILES,
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        name="aspirin",
        properties=CompoundProperties(
            molecular_weight=180.16,
            xlogp=1.2,
            tpsa=63.6,
            hbd_count=1,
            hba_count=4,
            rotatable_bonds=3,
            heavy_atom_count=13,
        ),
    )


def _make_caffeine() -> Compound:
    return Compound(
        cid=2519,
        smiles=CAFFEINE_SMILES,
        name="caffeine",
        properties=CompoundProperties(
            molecular_weight=194.19,
            xlogp=-0.07,
            tpsa=58.4,
            hbd_count=0,
            hba_count=6,
            rotatable_bonds=0,
            heavy_atom_count=14,
        ),
    )


def _make_cyclosporine() -> Compound:
    return Compound(
        cid=5284373,
        smiles=CYCLOSPORINE_SMILES,
        name="cyclosporine",
        properties=CompoundProperties(
            molecular_weight=1202.61,
            xlogp=2.9,
            tpsa=278.8,
            hbd_count=5,
            hba_count=23,
            rotatable_bonds=15,
            heavy_atom_count=88,
        ),
    )


def _make_no_smiles() -> Compound:
    return Compound(
        name="no_smiles_compound",
        properties=CompoundProperties(molecular_weight=180.0),
    )


# ---------------------------------------------------------------------------
# Compound.compute_descriptors()
# ---------------------------------------------------------------------------


class TestCompoundComputeDescriptors:
    def test_aspirin_returns_dict(self):
        compound = _make_aspirin()
        result = compound.compute_descriptors()
        assert isinstance(result, dict)

    def test_aspirin_200_plus_descriptors(self):
        compound = _make_aspirin()
        result = compound.compute_descriptors()
        assert len(result) >= 200

    def test_descriptors_stored_on_compound(self):
        compound = _make_aspirin()
        compound.compute_descriptors()
        assert len(compound.descriptors) >= 200

    def test_returns_same_as_stored(self):
        compound = _make_aspirin()
        result = compound.compute_descriptors()
        assert result is compound.descriptors

    def test_no_smiles_raises_value_error(self):
        compound = _make_no_smiles()
        with pytest.raises(ValueError, match="no SMILES"):
            compound.compute_descriptors()

    def test_values_are_numeric(self):
        compound = _make_aspirin()
        result = compound.compute_descriptors()
        for name, value in result.items():
            assert isinstance(value, (int, float)), f"{name} is not numeric"

    def test_caffeine_also_works(self):
        compound = _make_caffeine()
        result = compound.compute_descriptors()
        assert len(result) >= 200


# ---------------------------------------------------------------------------
# Compound.compute_fingerprints()
# ---------------------------------------------------------------------------


class TestCompoundComputeFingerprints:
    def test_default_types_morgan_and_maccs(self):
        compound = _make_aspirin()
        result = compound.compute_fingerprints()
        assert "morgan" in result
        assert "maccs" in result

    def test_fingerprints_stored_on_compound(self):
        compound = _make_aspirin()
        compound.compute_fingerprints()
        assert "morgan" in compound.fingerprints
        assert "maccs" in compound.fingerprints

    def test_returns_same_as_stored(self):
        compound = _make_aspirin()
        result = compound.compute_fingerprints()
        assert result is compound.fingerprints

    def test_custom_types(self):
        compound = _make_aspirin()
        result = compound.compute_fingerprints(types=["rdkit"])
        assert "rdkit" in result
        assert "morgan" not in result

    def test_morgan_result_has_required_keys(self):
        compound = _make_aspirin()
        result = compound.compute_fingerprints(types=["morgan"])
        fp = result["morgan"]
        assert "type" in fp
        assert "bits" in fp
        assert "num_bits" in fp
        assert "num_on_bits" in fp

    def test_custom_radius_and_n_bits(self):
        compound = _make_aspirin()
        result = compound.compute_fingerprints(types=["morgan"], radius=3, n_bits=1024)
        fp = result["morgan"]
        assert fp["num_bits"] == 1024

    def test_no_smiles_raises_value_error(self):
        compound = _make_no_smiles()
        with pytest.raises(ValueError, match="no SMILES"):
            compound.compute_fingerprints()

    def test_all_five_types(self):
        compound = _make_aspirin()
        fp_types = ["morgan", "maccs", "rdkit", "topological_torsion", "atom_pair"]
        result = compound.compute_fingerprints(types=fp_types)
        for fp_type in fp_types:
            assert fp_type in result


# ---------------------------------------------------------------------------
# Compound.check_drug_likeness()
# ---------------------------------------------------------------------------


class TestCompoundCheckDrugLikeness:
    def test_aspirin_returns_drug_likeness(self):
        from chemfuse.models.prediction import DrugLikeness

        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        assert isinstance(result, DrugLikeness)

    def test_druglikeness_stored_on_compound(self):
        from chemfuse.models.prediction import DrugLikeness

        compound = _make_aspirin()
        compound.check_drug_likeness()
        assert isinstance(compound.druglikeness, DrugLikeness)

    def test_returns_same_as_stored(self):
        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        assert result is compound.druglikeness

    def test_aspirin_lipinski_passes(self):
        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        assert result.lipinski.pass_filter is True

    def test_aspirin_veber_passes(self):
        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        assert result.veber.pass_filter is True

    def test_aspirin_overall_pass(self):
        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        # Aspirin should pass all five filters
        assert result.overall_pass is True

    def test_cyclosporine_fails_lipinski(self):
        compound = _make_cyclosporine()
        result = compound.check_drug_likeness()
        assert result.lipinski.pass_filter is False

    def test_cyclosporine_overall_fail(self):
        compound = _make_cyclosporine()
        result = compound.check_drug_likeness()
        assert result.overall_pass is False

    def test_all_five_filters_present(self):
        from chemfuse.models.prediction import FilterResult

        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        assert isinstance(result.lipinski, FilterResult)
        assert isinstance(result.veber, FilterResult)
        assert isinstance(result.ghose, FilterResult)
        assert isinstance(result.egan, FilterResult)
        assert isinstance(result.muegge, FilterResult)

    def test_summary_returns_dict_of_bools(self):
        compound = _make_aspirin()
        result = compound.check_drug_likeness()
        summary = result.summary()
        assert isinstance(summary, dict)
        core_keys = {"lipinski", "veber", "ghose", "egan", "muegge"}
        assert core_keys.issubset(set(summary.keys()))
        for val in summary.values():
            assert isinstance(val, bool)

    def test_uses_preloaded_properties(self):
        """Compound.check_drug_likeness() should use pre-populated properties."""
        compound = _make_aspirin()
        # Properties are set from PubChem data; the result should be deterministic
        result1 = compound.check_drug_likeness()
        result2 = compound.check_drug_likeness()
        assert result1.lipinski.pass_filter == result2.lipinski.pass_filter

    def test_no_smiles_still_works_with_properties(self):
        """Without SMILES, drug-likeness can still be evaluated from properties."""
        compound = Compound(
            name="no_smiles",
            smiles="",
            properties=CompoundProperties(
                molecular_weight=180.16,
                xlogp=1.2,
                tpsa=63.6,
                hbd_count=1,
                hba_count=4,
                rotatable_bonds=3,
                heavy_atom_count=13,
            ),
        )
        result = compound.check_drug_likeness()
        # Should be able to evaluate Lipinski and Veber from properties alone
        assert result.lipinski.pass_filter is True


# ---------------------------------------------------------------------------
# CompoundCollection.compute_all()
# ---------------------------------------------------------------------------


class TestCompoundCollectionComputeAll:
    def test_computes_descriptors_for_all(self):
        aspirin = _make_aspirin()
        caffeine = _make_caffeine()
        collection = CompoundCollection(compounds=[aspirin, caffeine])
        collection.compute_all(descriptors=True, fingerprints=False, druglikeness=False)
        assert len(aspirin.descriptors) >= 200
        assert len(caffeine.descriptors) >= 200

    def test_computes_druglikeness_for_all(self):
        from chemfuse.models.prediction import DrugLikeness

        aspirin = _make_aspirin()
        caffeine = _make_caffeine()
        collection = CompoundCollection(compounds=[aspirin, caffeine])
        collection.compute_all(descriptors=False, fingerprints=False, druglikeness=True)
        assert isinstance(aspirin.druglikeness, DrugLikeness)
        assert isinstance(caffeine.druglikeness, DrugLikeness)

    def test_computes_fingerprints_for_all(self):
        aspirin = _make_aspirin()
        caffeine = _make_caffeine()
        collection = CompoundCollection(compounds=[aspirin, caffeine])
        collection.compute_all(descriptors=False, fingerprints=True, druglikeness=False)
        assert "morgan" in aspirin.fingerprints
        assert "morgan" in caffeine.fingerprints

    def test_skips_no_smiles_gracefully(self):
        """Compounds without SMILES should not abort the whole batch."""
        aspirin = _make_aspirin()
        no_smiles = _make_no_smiles()
        collection = CompoundCollection(compounds=[aspirin, no_smiles])
        # Should not raise
        collection.compute_all(descriptors=True, fingerprints=False, druglikeness=False)
        # Aspirin descriptors computed; no_smiles skipped silently
        assert len(aspirin.descriptors) >= 200
        assert len(no_smiles.descriptors) == 0

    def test_progress_callback_called(self):
        aspirin = _make_aspirin()
        caffeine = _make_caffeine()
        collection = CompoundCollection(compounds=[aspirin, caffeine])
        calls: list[tuple[int, int, str]] = []
        collection.compute_all(
            descriptors=False,
            fingerprints=False,
            druglikeness=True,
            progress_callback=lambda done, total, name: calls.append((done, total, name)),
        )
        assert len(calls) == 2
        assert calls[0][0] == 1
        assert calls[1][0] == 2
        assert calls[0][1] == 2

    def test_empty_collection_no_error(self):
        collection = CompoundCollection()
        # Should complete without error
        collection.compute_all()

    def test_all_false_is_noop(self):
        aspirin = _make_aspirin()
        collection = CompoundCollection(compounds=[aspirin])
        collection.compute_all(descriptors=False, fingerprints=False, druglikeness=False)
        # Nothing computed
        assert len(aspirin.descriptors) == 0
        assert len(aspirin.fingerprints) == 0
        assert aspirin.druglikeness is None


# ---------------------------------------------------------------------------
# CompoundCollection.filter_by_druglikeness()
# ---------------------------------------------------------------------------


class TestCompoundCollectionFilterByDruglikeness:
    def _make_collection(self) -> CompoundCollection:
        aspirin = _make_aspirin()
        caffeine = _make_caffeine()
        cyclosporine = _make_cyclosporine()
        return CompoundCollection(compounds=[aspirin, caffeine, cyclosporine])

    def test_lipinski_filter_removes_cyclosporine(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(lipinski=True)
        names = [c.name for c in result]
        assert "aspirin" in names
        assert "caffeine" in names
        assert "cyclosporine" not in names

    def test_returns_compound_collection(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(lipinski=True)
        assert isinstance(result, CompoundCollection)

    def test_no_filters_keeps_all(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness()
        assert len(result) == 3

    def test_multiple_filters_applied_together(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(lipinski=True, veber=True)
        # Cyclosporine fails both Lipinski and Veber
        names = [c.name for c in result]
        assert "cyclosporine" not in names

    def test_pre_computed_druglikeness_is_reused(self):
        """If druglikeness already set, filter should use it without recomputing."""

        aspirin = _make_aspirin()
        # Pre-compute and mark drug-likeness
        aspirin.check_drug_likeness()
        assert aspirin.druglikeness is not None

        stored_dl = aspirin.druglikeness
        collection = CompoundCollection(compounds=[aspirin])
        collection.filter_by_druglikeness(lipinski=True)
        # Object should be the same (not recomputed)
        assert aspirin.druglikeness is stored_dl

    def test_empty_collection_returns_empty(self):
        collection = CompoundCollection()
        result = collection.filter_by_druglikeness(lipinski=True)
        assert len(result) == 0

    def test_veber_filter(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(veber=True)
        names = [c.name for c in result]
        # Cyclosporine has TPSA=278.8, far above Veber's 140 limit
        assert "cyclosporine" not in names
        assert "aspirin" in names

    def test_ghose_filter(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(ghose=True)
        names = [c.name for c in result]
        # Cyclosporine MW=1202.61 far outside 160-480 range
        assert "cyclosporine" not in names

    def test_egan_filter(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(egan=True)
        names = [c.name for c in result]
        # Cyclosporine TPSA=278.8 > 131.6
        assert "cyclosporine" not in names
        assert "aspirin" in names

    def test_muegge_filter(self):
        collection = self._make_collection()
        result = collection.filter_by_druglikeness(muegge=True)
        names = [c.name for c in result]
        # Cyclosporine MW=1202.61 > 600
        assert "cyclosporine" not in names
