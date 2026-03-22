"""Tests for analyze/diversity.py (CF-E06: MaxMin Diversity Picking)."""

from __future__ import annotations

import pytest

from chemfuse.core.exceptions import OptionalDependencyError

# ---------------------------------------------------------------------------
# Test SMILES sets
# ---------------------------------------------------------------------------

# Similar aliphatic alcohols (low diversity)
SIMILAR_SMILES = [
    "CCO",        # ethanol
    "CCCO",       # propanol
    "CCCCO",      # butanol
    "CCCCCO",     # pentanol
    "CCCCCCO",    # hexanol
]

# Diverse compounds across different scaffolds / pharmacophores
DIVERSE_SMILES = [
    "c1ccccc1",        # benzene
    "C1CCCCC1",        # cyclohexane
    "CC(=O)O",         # acetic acid
    "c1ccncc1",        # pyridine
    "C1CCNCC1",        # piperidine
]

# Mixed: 5 similar + 15 diverse (used for diversity comparison tests)
MIXED_SMILES = SIMILAR_SMILES + [
    "c1ccccc1",
    "C1CCCCC1",
    "CC(=O)O",
    "c1ccncc1",
    "C1CCNCC1",
    "c1ccc(O)cc1",     # phenol
    "CC(N)=O",         # acetamide
    "c1ccsc1",         # thiophene
    "C1CCOC1",         # THF
    "c1ccc2ncccc2c1",  # quinoline
    "CC(=O)Nc1ccccc1", # acetanilide
    "c1ccc(Cl)cc1",    # chlorobenzene
    "CCN(CC)CC",       # triethylamine
    "OC(=O)c1ccccc1",  # benzoic acid
    "C1CCCCC1N",       # cyclohexylamine
]

# Twenty structurally diverse compounds
DIVERSE_20 = [
    "c1ccccc1",
    "C1CCCCC1",
    "CC(=O)O",
    "c1ccncc1",
    "C1CCNCC1",
    "c1ccc(O)cc1",
    "CC(N)=O",
    "c1ccsc1",
    "C1CCOC1",
    "c1ccc2ncccc2c1",
    "CC(=O)Nc1ccccc1",
    "c1ccc(Cl)cc1",
    "CCN(CC)CC",
    "OC(=O)c1ccccc1",
    "C1CCCCC1N",
    "CC(=O)c1ccccc1",
    "c1ccnc(N)c1",
    "C(F)(F)(F)c1ccccc1",
    "c1cnc2ccccc2n1",
    "CCO",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_compound(smiles: str):
    """Create a minimal Compound-like object for CompoundCollection tests."""
    from chemfuse.models.compound import Compound
    return Compound(smiles=smiles, name=smiles)


def _make_collection(smiles_list: list[str]):
    """Build a CompoundCollection from a list of SMILES strings."""
    from chemfuse.models.collection import CompoundCollection
    compounds = [_make_compound(s) for s in smiles_list]
    return CompoundCollection(compounds=compounds)


# ---------------------------------------------------------------------------
# Tests for maxmin_pick
# ---------------------------------------------------------------------------

class TestMaxminPick:
    def test_returns_n_pick_distinct_indices(self):
        """maxmin_pick returns exactly n_pick distinct indices from a 20-compound set."""
        from chemfuse.analyze.diversity import maxmin_pick

        n_pick = 5
        picks = maxmin_pick(DIVERSE_20, n_pick=n_pick)

        assert len(picks) == n_pick
        assert len(set(picks)) == n_pick  # all distinct
        assert all(0 <= i < len(DIVERSE_20) for i in picks)

    def test_picked_subset_more_diverse_than_similar_subset(self):
        """MaxMin picks from mixed set are more diverse than a random similar subset."""
        from chemfuse.analyze.diversity import diversity_score, maxmin_pick

        n_pick = 5
        picks = maxmin_pick(MIXED_SMILES, n_pick=n_pick)
        picked_smiles = [MIXED_SMILES[i] for i in picks]

        # Random subset of 5 similar aliphatic alcohols
        similar_subset = SIMILAR_SMILES[:5]

        picked_diversity = diversity_score(picked_smiles)
        similar_diversity = diversity_score(similar_subset)

        assert picked_diversity > similar_diversity, (
            f"Expected picked diversity ({picked_diversity:.3f}) > "
            f"similar diversity ({similar_diversity:.3f})"
        )

    def test_n_pick_greater_than_len_returns_all_valid_indices(self):
        """When n_pick >= len(smiles_list), all valid indices are returned."""
        from chemfuse.analyze.diversity import maxmin_pick

        smiles = ["CCO", "CCCO", "c1ccccc1"]
        picks = maxmin_pick(smiles, n_pick=10)

        assert sorted(picks) == [0, 1, 2]

    def test_n_pick_equals_len_returns_all(self):
        """When n_pick == len(smiles_list), all valid indices are returned."""
        from chemfuse.analyze.diversity import maxmin_pick

        smiles = ["CCO", "CCCO", "c1ccccc1"]
        picks = maxmin_pick(smiles, n_pick=3)

        assert sorted(picks) == [0, 1, 2]

    def test_n_pick_one_returns_single_index(self):
        """n_pick=1 returns exactly one valid index."""
        from chemfuse.analyze.diversity import maxmin_pick

        picks = maxmin_pick(DIVERSE_20, n_pick=1)

        assert len(picks) == 1
        assert 0 <= picks[0] < len(DIVERSE_20)

    def test_empty_list_returns_empty(self):
        """Empty SMILES list returns empty picks."""
        from chemfuse.analyze.diversity import maxmin_pick

        picks = maxmin_pick([], n_pick=5)
        assert picks == []

    def test_invalid_smiles_excluded(self):
        """Invalid SMILES are silently excluded; valid ones are returned."""
        from chemfuse.analyze.diversity import maxmin_pick

        smiles = ["CCO", "INVALID_SMILES###", "c1ccccc1", "CCCO"]
        picks = maxmin_pick(smiles, n_pick=3)

        # Only 3 valid SMILES available, so all valid indices returned
        valid_indices = [0, 2, 3]
        assert sorted(picks) == valid_indices

    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is unavailable."""
        from unittest.mock import patch

        with patch("chemfuse.analyze.diversity._RDKIT_AVAILABLE", False):
            from chemfuse.analyze import diversity
            with pytest.raises(OptionalDependencyError):
                diversity.maxmin_pick(["CCO", "c1ccccc1"], n_pick=1)

    def test_maccs_fingerprint_type(self):
        """maxmin_pick works with maccs fingerprint type."""
        from chemfuse.analyze.diversity import maxmin_pick

        picks = maxmin_pick(DIVERSE_20, n_pick=5, fp_type="maccs")
        assert len(picks) == 5
        assert len(set(picks)) == 5

    def test_rdkit_fingerprint_type(self):
        """maxmin_pick works with rdkit fingerprint type."""
        from chemfuse.analyze.diversity import maxmin_pick

        picks = maxmin_pick(DIVERSE_20, n_pick=5, fp_type="rdkit")
        assert len(picks) == 5
        assert len(set(picks)) == 5


# ---------------------------------------------------------------------------
# Tests for diversity_score
# ---------------------------------------------------------------------------

class TestDiversityScore:
    def test_identical_smiles_returns_zero(self):
        """All identical SMILES yield diversity_score = 0.0."""
        from chemfuse.analyze.diversity import diversity_score

        smiles = ["CCO", "CCO", "CCO", "CCO"]
        score = diversity_score(smiles)

        assert score == pytest.approx(0.0, abs=1e-6)

    def test_very_different_smiles_returns_high_score(self):
        """Structurally diverse compounds yield diversity_score > 0.5."""
        from chemfuse.analyze.diversity import diversity_score

        score = diversity_score(DIVERSE_SMILES)

        assert score > 0.5, f"Expected diversity > 0.5, got {score:.3f}"

    def test_diverse_greater_than_similar(self):
        """Diverse compounds score higher than similar aliphatic alcohols."""
        from chemfuse.analyze.diversity import diversity_score

        diverse = diversity_score(DIVERSE_SMILES)
        similar = diversity_score(SIMILAR_SMILES)

        assert diverse > similar

    def test_single_compound_returns_zero(self):
        """A single-compound list returns 0.0 (no pairwise distances)."""
        from chemfuse.analyze.diversity import diversity_score

        score = diversity_score(["CCO"])
        assert score == pytest.approx(0.0, abs=1e-6)

    def test_empty_list_returns_zero(self):
        """Empty list returns 0.0."""
        from chemfuse.analyze.diversity import diversity_score

        score = diversity_score([])
        assert score == pytest.approx(0.0, abs=1e-6)

    def test_returns_float_in_unit_interval(self):
        """diversity_score always returns a float in [0.0, 1.0]."""
        from chemfuse.analyze.diversity import diversity_score

        score = diversity_score(DIVERSE_20)
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is unavailable."""
        from unittest.mock import patch

        with patch("chemfuse.analyze.diversity._RDKIT_AVAILABLE", False):
            from chemfuse.analyze import diversity
            with pytest.raises(OptionalDependencyError):
                diversity.diversity_score(["CCO", "c1ccccc1"])


# ---------------------------------------------------------------------------
# Tests for CompoundCollection.pick_diverse
# ---------------------------------------------------------------------------

class TestCompoundCollectionPickDiverse:
    def test_returns_compound_collection_of_correct_size(self):
        """pick_diverse returns a CompoundCollection with exactly n compounds."""
        collection = _make_collection(DIVERSE_20)
        picked = collection.pick_diverse(n=5)

        assert len(picked) == 5

    def test_pick_diverse_n_greater_than_size_returns_all(self):
        """pick_diverse with n > collection size returns all compounds."""
        collection = _make_collection(DIVERSE_SMILES)
        picked = collection.pick_diverse(n=100)

        assert len(picked) == len(DIVERSE_SMILES)

    def test_returns_compound_collection_instance(self):
        """pick_diverse returns a CompoundCollection, not a plain list."""
        from chemfuse.models.collection import CompoundCollection

        collection = _make_collection(DIVERSE_20)
        picked = collection.pick_diverse(n=5)

        assert isinstance(picked, CompoundCollection)

    def test_picked_compounds_are_from_original(self):
        """All picked compounds exist in the original collection."""
        collection = _make_collection(DIVERSE_20)
        picked = collection.pick_diverse(n=5)

        original_smiles = {c.smiles for c in collection}
        picked_smiles = {c.smiles for c in picked}
        assert picked_smiles.issubset(original_smiles)

    def test_empty_collection_returns_empty(self):
        """pick_diverse on empty collection returns empty collection."""
        collection = _make_collection([])
        picked = collection.pick_diverse(n=5)

        assert len(picked) == 0

    def test_pick_diverse_one_returns_single_compound(self):
        """pick_diverse(n=1) returns exactly one compound."""
        collection = _make_collection(DIVERSE_20)
        picked = collection.pick_diverse(n=1)

        assert len(picked) == 1


# ---------------------------------------------------------------------------
# Tests for CompoundCollection.diversity_score
# ---------------------------------------------------------------------------

class TestCompoundCollectionDiversityScore:
    def test_returns_float_in_unit_interval(self):
        """diversity_score returns a float in [0.0, 1.0]."""
        collection = _make_collection(DIVERSE_20)
        score = collection.diversity_score()

        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_identical_compounds_score_zero(self):
        """Collection of identical compounds scores 0.0."""
        collection = _make_collection(["CCO", "CCO", "CCO"])
        score = collection.diversity_score()

        assert score == pytest.approx(0.0, abs=1e-6)

    def test_diverse_compounds_score_above_zero(self):
        """Diverse collection scores > 0.0."""
        collection = _make_collection(DIVERSE_SMILES)
        score = collection.diversity_score()

        assert score > 0.0

    def test_empty_collection_returns_zero(self):
        """Empty collection returns 0.0."""
        collection = _make_collection([])
        score = collection.diversity_score()

        assert score == pytest.approx(0.0, abs=1e-6)

    def test_single_compound_returns_zero(self):
        """Single-compound collection returns 0.0."""
        collection = _make_collection(["CCO"])
        score = collection.diversity_score()

        assert score == pytest.approx(0.0, abs=1e-6)

    def test_diverse_greater_than_similar(self):
        """Diverse collection scores higher than similar aliphatic alcohols."""
        diverse_coll = _make_collection(DIVERSE_SMILES)
        similar_coll = _make_collection(SIMILAR_SMILES)

        assert diverse_coll.diversity_score() > similar_coll.diversity_score()
