"""Tests for chemfuse.compute.fingerprints."""

from __future__ import annotations

import pytest

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(c1=O)C"
INVALID_SMILES = "not_a_smiles_string!!!"


class TestComputeFingerprint:
    """Tests for compute_fingerprint()."""

    @pytest.mark.parametrize("fp_type", [
        "morgan", "maccs", "rdkit", "topological_torsion", "atom_pair"
    ])
    def test_all_fingerprint_types_return_dict(self, fp_type):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type=fp_type)
        assert result is not None
        assert isinstance(result, dict)

    @pytest.mark.parametrize("fp_type", [
        "morgan", "maccs", "rdkit", "topological_torsion", "atom_pair"
    ])
    def test_result_has_required_keys(self, fp_type):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type=fp_type)
        assert result is not None
        assert "type" in result
        assert "bits" in result
        assert "num_bits" in result
        assert "num_on_bits" in result

    def test_morgan_default_2048_bits(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type="morgan")
        assert result is not None
        assert result["num_bits"] == 2048

    def test_maccs_166_bits(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type="maccs")
        assert result is not None
        assert result["num_bits"] == 167  # MACCS is 0-indexed, so 167 positions

    def test_custom_radius(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result_r2 = compute_fingerprint(ASPIRIN_SMILES, fp_type="morgan", radius=2)
        result_r3 = compute_fingerprint(ASPIRIN_SMILES, fp_type="morgan", radius=3)
        # Different radii should generally produce different fingerprints
        assert result_r2 is not None
        assert result_r3 is not None
        # They may differ (different chemical environment captured)

    def test_custom_n_bits(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type="morgan", n_bits=1024)
        assert result is not None
        assert result["num_bits"] == 1024

    def test_bits_is_list_of_ints(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type="morgan")
        assert result is not None
        assert isinstance(result["bits"], list)
        for bit in result["bits"]:
            assert isinstance(bit, int)

    def test_invalid_smiles_returns_none(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(INVALID_SMILES, fp_type="morgan")
        assert result is None

    def test_invalid_smiles_logs_warning(self, caplog):
        import logging

        from chemfuse.compute.fingerprints import compute_fingerprint
        with caplog.at_level(logging.WARNING, logger="chemfuse.compute.fingerprints"):
            compute_fingerprint(INVALID_SMILES)
        assert any("Invalid SMILES" in r.message for r in caplog.records)

    def test_unknown_fp_type_raises_value_error(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        with pytest.raises(ValueError, match="Unknown fingerprint type"):
            compute_fingerprint(ASPIRIN_SMILES, fp_type="nonexistent_type")

    def test_rdkit_not_installed_raises_import_error(self):
        import chemfuse.compute.fingerprints as fp_mod
        original = fp_mod.RDKIT_AVAILABLE
        fp_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                fp_mod.compute_fingerprint(ASPIRIN_SMILES)
        finally:
            fp_mod.RDKIT_AVAILABLE = original

    def test_type_field_matches_requested(self):
        from chemfuse.compute.fingerprints import compute_fingerprint
        result = compute_fingerprint(ASPIRIN_SMILES, fp_type="maccs")
        assert result is not None
        assert result["type"] == "maccs"


class TestTanimotoSimilarity:
    """Tests for tanimoto_similarity()."""

    def test_identical_molecules_similarity_one(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        sim = tanimoto_similarity(ASPIRIN_SMILES, ASPIRIN_SMILES)
        assert sim == pytest.approx(1.0)

    def test_different_molecules_similarity_less_than_one(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        sim = tanimoto_similarity(ASPIRIN_SMILES, CAFFEINE_SMILES)
        assert 0.0 <= sim < 1.0

    def test_similarity_is_float(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        sim = tanimoto_similarity(ASPIRIN_SMILES, CAFFEINE_SMILES)
        assert isinstance(sim, float)

    def test_similarity_in_range(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        sim = tanimoto_similarity(ASPIRIN_SMILES, CAFFEINE_SMILES)
        assert 0.0 <= sim <= 1.0

    @pytest.mark.parametrize("fp_type", ["morgan", "maccs", "rdkit"])
    def test_similarity_all_fp_types(self, fp_type):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        sim = tanimoto_similarity(ASPIRIN_SMILES, CAFFEINE_SMILES, fp_type=fp_type)
        assert 0.0 <= sim <= 1.0

    def test_invalid_smiles_a_raises(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        with pytest.raises(ValueError):
            tanimoto_similarity(INVALID_SMILES, ASPIRIN_SMILES)

    def test_invalid_smiles_b_raises(self):
        from chemfuse.compute.fingerprints import tanimoto_similarity
        with pytest.raises(ValueError):
            tanimoto_similarity(ASPIRIN_SMILES, INVALID_SMILES)

    def test_rdkit_not_installed_raises(self):
        import chemfuse.compute.fingerprints as fp_mod
        original = fp_mod.RDKIT_AVAILABLE
        fp_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                fp_mod.tanimoto_similarity(ASPIRIN_SMILES, CAFFEINE_SMILES)
        finally:
            fp_mod.RDKIT_AVAILABLE = original


class TestBulkTanimoto:
    """Tests for bulk_tanimoto()."""

    def test_returns_list_of_tuples(self):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        targets = [CAFFEINE_SMILES, "CC(=O)O"]
        results = bulk_tanimoto(ASPIRIN_SMILES, targets)
        assert isinstance(results, list)
        for item in results:
            assert isinstance(item, tuple)
            assert len(item) == 2
            smiles, sim = item
            assert isinstance(smiles, str)
            assert isinstance(sim, float)

    def test_results_sorted_descending(self):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        targets = [CAFFEINE_SMILES, "CC(=O)O", "c1ccccc1"]
        results = bulk_tanimoto(ASPIRIN_SMILES, targets)
        similarities = [r[1] for r in results]
        assert similarities == sorted(similarities, reverse=True)

    def test_invalid_target_skipped(self):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        targets = [CAFFEINE_SMILES, INVALID_SMILES, "CC(=O)O"]
        results = bulk_tanimoto(ASPIRIN_SMILES, targets)
        # Invalid SMILES should be silently skipped
        assert len(results) == 2

    def test_invalid_query_raises(self):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        with pytest.raises(ValueError, match="Invalid query"):
            bulk_tanimoto(INVALID_SMILES, [ASPIRIN_SMILES])

    def test_empty_target_list(self):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        results = bulk_tanimoto(ASPIRIN_SMILES, [])
        assert results == []

    @pytest.mark.parametrize("fp_type", ["morgan", "maccs", "rdkit", "topological_torsion", "atom_pair"])
    def test_all_fp_types(self, fp_type):
        from chemfuse.compute.fingerprints import bulk_tanimoto
        results = bulk_tanimoto(ASPIRIN_SMILES, [CAFFEINE_SMILES], fp_type=fp_type)
        assert len(results) == 1
        assert 0.0 <= results[0][1] <= 1.0

    def test_rdkit_not_installed_raises(self):
        import chemfuse.compute.fingerprints as fp_mod
        original = fp_mod.RDKIT_AVAILABLE
        fp_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                fp_mod.bulk_tanimoto(ASPIRIN_SMILES, [CAFFEINE_SMILES])
        finally:
            fp_mod.RDKIT_AVAILABLE = original
