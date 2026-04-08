"""Tests for analyze/similarity.py (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.core.exceptions import OptionalDependencyError

SMILES_A = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
SMILES_B = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"  # caffeine


def _make_fp_mock(bits=2048):
    fp = MagicMock()
    return fp


def _patch_rdkit(available=True):
    """Return a patch context for rdkit availability in the similarity module."""
    return patch("chemfuse.analyze.similarity._RDKIT_AVAILABLE", available)


class TestTanimotoSimilarity:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.tanimoto_similarity(SMILES_A, SMILES_B)

    def test_identical_smiles_returns_one(self):
        """tanimoto_similarity for identical SMILES returns 1.0."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_morgan_gen.GetFingerprint.return_value = fp
            mock_ds.TanimotoSimilarity.return_value = 1.0

            from chemfuse.analyze import similarity
            result = similarity.tanimoto_similarity(SMILES_A, SMILES_A)
            assert result == 1.0

    def test_different_smiles_returns_float(self):
        """tanimoto_similarity returns a float in [0, 1]."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_morgan_gen.GetFingerprint.return_value = fp
            mock_ds.TanimotoSimilarity.return_value = 0.35

            from chemfuse.analyze import similarity
            result = similarity.tanimoto_similarity(SMILES_A, SMILES_B)
            assert isinstance(result, float)
            assert 0.0 <= result <= 1.0

    def test_invalid_smiles_raises_value_error(self):
        """tanimoto_similarity raises ValueError for invalid SMILES."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:

            mock_chem.MolFromSmiles.return_value = None

            from chemfuse.analyze import similarity
            with pytest.raises(ValueError, match="Invalid SMILES"):
                similarity.tanimoto_similarity("BADSMILES", SMILES_B)


class TestBulkTanimoto:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.bulk_tanimoto(SMILES_A, [SMILES_B])

    def test_returns_sorted_descending(self):
        """bulk_tanimoto returns results sorted by descending similarity."""
        sims = [0.5, 0.9, 0.1]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_morgan_gen.GetFingerprint.return_value = fp
            mock_ds.TanimotoSimilarity.side_effect = sims

            from chemfuse.analyze import similarity
            results = similarity.bulk_tanimoto(SMILES_A, [SMILES_A, SMILES_B, "CCO"])

        assert len(results) == 3
        scores = [r[1] for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_empty_targets_returns_empty(self):
        """bulk_tanimoto with empty targets returns empty list."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_morgan_gen.GetFingerprint.return_value = _make_fp_mock()

            from chemfuse.analyze import similarity
            results = similarity.bulk_tanimoto(SMILES_A, [])

        assert results == []


class TestTanimotoMatrix:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.tanimoto_matrix([SMILES_A])

    def test_matrix_shape(self):
        """tanimoto_matrix returns n x n array."""
        smiles_list = [SMILES_A, SMILES_B, "CCO"]
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_morgan_gen.GetFingerprint.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.4

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        assert isinstance(mat, np.ndarray)
        assert mat.shape == (3, 3)

    def test_matrix_diagonal_is_one(self):
        """Diagonal of tanimoto_matrix should be 1.0 (self-similarity)."""
        smiles_list = [SMILES_A, SMILES_B]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_morgan_gen.GetFingerprint.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.35

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        assert mat[0, 0] == pytest.approx(1.0)
        assert mat[1, 1] == pytest.approx(1.0)

    def test_matrix_is_symmetric(self):
        """tanimoto_matrix must be symmetric."""
        smiles_list = [SMILES_A, SMILES_B, "CCO"]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_morgan_gen.GetFingerprint.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.4

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        np.testing.assert_array_equal(mat, mat.T)


class TestFindNearestNeighbors:
    def test_returns_top_n(self):
        """find_nearest_neighbors returns at most n results."""
        targets = [SMILES_A, SMILES_B, "CCO", "c1ccccc1", "CCC"]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity._MORGAN_GEN") as mock_morgan_gen:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_morgan_gen.GetFingerprint.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.5

            from chemfuse.analyze import similarity
            results = similarity.find_nearest_neighbors(SMILES_A, targets, n=3)

        assert len(results) == 3


# ---------------------------------------------------------------------------
# CF-E04: Substructure Search (SMARTS)
# ---------------------------------------------------------------------------

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ETHANOL = "CCO"
SULFAMETHOXAZOLE = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
NICOTINE = "CN1CCCC1c1cccnc1"
BENZENE_SMARTS = "c1ccccc1"
SULFONAMIDE_SMARTS = "S(=O)(=O)N"
PYRIDINE_SMARTS = "c1ccncc1"


def _make_mol_mock(has_match: bool = True, match_atoms: tuple = (0, 1, 2)):
    """Create a mock RDKit Mol that returns a fixed substructure match result."""
    mol = MagicMock()
    mol.HasSubstructMatch.return_value = has_match
    mol.GetSubstructMatches.return_value = (match_atoms,) if has_match else ()
    return mol


class TestSubstructureSearch:
    """Tests for substructure_search (CF-E04)."""

    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.substructure_search([ASPIRIN], BENZENE_SMARTS)

    def test_invalid_smarts_raises_value_error(self):
        """Invalid SMARTS pattern raises ValueError."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = None

            from chemfuse.analyze import similarity
            with pytest.raises(ValueError, match="Invalid SMARTS"):
                similarity.substructure_search([ASPIRIN], "BADSMARTS!!!")

    def test_benzene_matches_aspirin_true(self):
        """Benzene ring SMARTS matches aspirin (True)."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            pattern = MagicMock()
            mock_chem.MolFromSmarts.return_value = pattern
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(has_match=True)

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([ASPIRIN], BENZENE_SMARTS)

        assert results == [True]

    def test_benzene_does_not_match_ethanol_false(self):
        """Benzene ring SMARTS does not match ethanol (False)."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            pattern = MagicMock()
            mock_chem.MolFromSmarts.return_value = pattern
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(has_match=False)

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([ETHANOL], BENZENE_SMARTS)

        assert results == [False]

    def test_sulfonamide_matches_sulfamethoxazole(self):
        """Sulfonamide SMARTS matches sulfamethoxazole (True)."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(has_match=True)

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([SULFAMETHOXAZOLE], SULFONAMIDE_SMARTS)

        assert results == [True]

    def test_pyridine_matches_nicotine(self):
        """Pyridine SMARTS matches nicotine (True)."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(has_match=True)

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([NICOTINE], PYRIDINE_SMARTS)

        assert results == [True]

    def test_invalid_smiles_returns_false_not_exception(self):
        """Invalid SMILES in the list returns False rather than raising."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            # First call (valid mol) returns a match, second (invalid) returns None
            mock_chem.MolFromSmiles.side_effect = [_make_mol_mock(True), None]

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([ASPIRIN, "NOTASMILES"], BENZENE_SMARTS)

        assert results == [True, False]

    def test_empty_smiles_list_returns_empty_list(self):
        """Empty smiles_list returns empty list."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()

            from chemfuse.analyze import similarity
            results = similarity.substructure_search([], BENZENE_SMARTS)

        assert results == []

    def test_multiple_compounds_mixed_results(self):
        """substructure_search returns correct True/False for a mixed list."""
        match_sequence = [True, False, True]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.side_effect = [
                _make_mol_mock(m) for m in match_sequence
            ]

            from chemfuse.analyze import similarity
            results = similarity.substructure_search(
                [ASPIRIN, ETHANOL, SULFAMETHOXAZOLE], BENZENE_SMARTS
            )

        assert results == match_sequence


class TestSubstructureMatchAtoms:
    """Tests for substructure_match_atoms (CF-E04)."""

    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.substructure_match_atoms(ASPIRIN, BENZENE_SMARTS)

    def test_invalid_smarts_raises_value_error(self):
        """Invalid SMARTS pattern raises ValueError."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = None

            from chemfuse.analyze import similarity
            with pytest.raises(ValueError, match="Invalid SMARTS"):
                similarity.substructure_match_atoms(ASPIRIN, "BADSMARTS!!!")

    def test_invalid_smiles_raises_value_error(self):
        """Invalid SMILES raises ValueError."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = None

            from chemfuse.analyze import similarity
            with pytest.raises(ValueError, match="Invalid SMILES"):
                similarity.substructure_match_atoms("BADSMILES", BENZENE_SMARTS)

    def test_returns_correct_atom_indices(self):
        """Returns list of tuples with atom indices for each match."""
        expected_matches = ((0, 1, 2, 3, 4, 5),)

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mol = MagicMock()
            mol.GetSubstructMatches.return_value = expected_matches
            mock_chem.MolFromSmiles.return_value = mol

            from chemfuse.analyze import similarity
            result = similarity.substructure_match_atoms(ASPIRIN, BENZENE_SMARTS)

        assert result == list(expected_matches)

    def test_no_match_returns_empty_list(self):
        """Returns empty list when pattern does not match the molecule."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mol = MagicMock()
            mol.GetSubstructMatches.return_value = ()
            mock_chem.MolFromSmiles.return_value = mol

            from chemfuse.analyze import similarity
            result = similarity.substructure_match_atoms(ETHANOL, BENZENE_SMARTS)

        assert result == []

    def test_multiple_matches_returns_all(self):
        """Returns all match tuples when pattern appears multiple times."""
        # Naphthalene has two fused rings; two benzene matches
        expected = ((0, 1, 2, 3, 4, 5), (4, 5, 6, 7, 8, 9))

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mol = MagicMock()
            mol.GetSubstructMatches.return_value = expected
            mock_chem.MolFromSmiles.return_value = mol

            from chemfuse.analyze import similarity
            result = similarity.substructure_match_atoms("c1ccc2ccccc2c1", BENZENE_SMARTS)

        assert len(result) == 2
        assert result[0] == (0, 1, 2, 3, 4, 5)
        assert result[1] == (4, 5, 6, 7, 8, 9)


class TestFilterBySubstructure:
    """Tests for CompoundCollection.filter_by_substructure (CF-E04)."""

    def _make_collection(self, smiles_list: list[str]) -> object:
        """Create a CompoundCollection from a list of SMILES strings."""
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound

        compounds = [Compound(smiles=smi, name=f"cpd_{i}") for i, smi in enumerate(smiles_list)]
        return CompoundCollection(compounds=compounds)

    def test_filter_returns_matching_compounds_only(self):
        """filter_by_substructure returns only compounds with True matches."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            # aspirin matches, ethanol does not
            mock_chem.MolFromSmiles.side_effect = [
                _make_mol_mock(True),
                _make_mol_mock(False),
            ]

            collection = self._make_collection([ASPIRIN, ETHANOL])
            result = collection.filter_by_substructure(BENZENE_SMARTS)

        assert len(result) == 1
        assert result.compounds[0].smiles == ASPIRIN

    def test_filter_no_matches_returns_empty_collection(self):
        """filter_by_substructure returns empty collection when no compounds match."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(False)

            collection = self._make_collection([ETHANOL, "CCO", "CCCO"])
            result = collection.filter_by_substructure(BENZENE_SMARTS)

        assert len(result) == 0

    def test_filter_all_match_returns_full_collection(self):
        """filter_by_substructure returns full collection when all compounds match."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(True)

            collection = self._make_collection([ASPIRIN, SULFAMETHOXAZOLE, NICOTINE])
            result = collection.filter_by_substructure(BENZENE_SMARTS)

        assert len(result) == 3

    def test_filter_preserves_metadata(self):
        """filter_by_substructure preserves query, sources, timestamp metadata."""
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound

        compounds = [Compound(smiles=ASPIRIN, name="aspirin")]
        collection = CompoundCollection(
            compounds=compounds,
            query="test_query",
            sources=["chembl"],
        )

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolFromSmiles.return_value = _make_mol_mock(True)

            result = collection.filter_by_substructure(BENZENE_SMARTS)

        assert result.query == "test_query"
        assert result.sources == ["chembl"]

    def test_invalid_smarts_propagates_value_error(self):
        """filter_by_substructure propagates ValueError for invalid SMARTS."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = None

            collection = self._make_collection([ASPIRIN])
            with pytest.raises(ValueError, match="Invalid SMARTS"):
                collection.filter_by_substructure("BADSMARTS!!!")

    def test_filter_empty_collection_returns_empty(self):
        """filter_by_substructure on an empty collection returns empty collection."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:
            mock_chem.MolFromSmarts.return_value = MagicMock()

            collection = self._make_collection([])
            result = collection.filter_by_substructure(BENZENE_SMARTS)

        assert len(result) == 0
