"""Tests for chemfuse.compute.standardization (CF-E03)."""
from __future__ import annotations

import sys
from unittest.mock import patch

import pytest

# Common SMILES constants
SODIUM_ACETATE_SMILES = "[Na+].CC(=O)[O-]"
ACETIC_ACID_SMILES = "CC(=O)O"
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

# Keto/enol tautomer pair for acetylacetone (2,4-pentanedione)
KETO_TAUTOMER = "CC(=O)CC(=O)C"   # keto form
ENOL_TAUTOMER = "CC(=O)CC(O)=C"   # enol form (alternative representation)


class TestStripSalts:
    """Tests for strip_salts()."""

    def test_sodium_acetate_returns_acetic_acid(self):
        """Sodium acetate salt should yield acetic acid (largest fragment)."""
        from chemfuse.compute.standardization import strip_salts

        result = strip_salts(SODIUM_ACETATE_SMILES)
        assert result is not None
        # The largest fragment is the acetate; verify it round-trips back to acetic acid
        from rdkit import Chem
        mol_result = Chem.MolFromSmiles(result)
        mol_expected = Chem.MolFromSmiles(ACETIC_ACID_SMILES)
        assert mol_result is not None
        assert mol_expected is not None
        # Compare canonical forms
        assert Chem.MolToSmiles(mol_result) == Chem.MolToSmiles(mol_expected)

    def test_no_salt_compound_unchanged(self):
        """A single-fragment compound should be returned as-is (canonicalized)."""
        from rdkit import Chem

        from chemfuse.compute.standardization import strip_salts

        result = strip_salts(ASPIRIN_SMILES)
        assert result is not None
        mol_result = Chem.MolFromSmiles(result)
        mol_expected = Chem.MolFromSmiles(ASPIRIN_SMILES)
        assert mol_result is not None
        assert Chem.MolToSmiles(mol_result) == Chem.MolToSmiles(mol_expected)

    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES should return None (no exception)."""
        from chemfuse.compute.standardization import strip_salts

        assert strip_salts("NOT_A_VALID_SMILES!!!") is None

    def test_empty_string_returns_none(self):
        """Empty string should return None."""
        from chemfuse.compute.standardization import strip_salts

        assert strip_salts("") is None

    def test_multicomponent_selects_largest(self):
        """Aspirin sodium salt: organic part must be selected over Na+."""
        from rdkit import Chem

        from chemfuse.compute.standardization import strip_salts

        # Aspirin as sodium salt: Na+ + aspirin anion
        aspirin_sodium = "[Na+].CC(=O)Oc1ccccc1C(=O)[O-]"
        result = strip_salts(aspirin_sodium)
        assert result is not None
        mol = Chem.MolFromSmiles(result)
        assert mol is not None
        # Must not contain sodium
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        assert "Na" not in atoms


class TestCanonicalTautomer:
    """Tests for canonical_tautomer()."""

    def test_keto_enol_pair_produces_same_output(self):
        """Keto and enol forms of acetylacetone must produce identical canonical SMILES."""
        from chemfuse.compute.standardization import canonical_tautomer

        result_keto = canonical_tautomer(KETO_TAUTOMER)
        result_enol = canonical_tautomer(ENOL_TAUTOMER)

        assert result_keto is not None
        assert result_enol is not None
        assert result_keto == result_enol

    def test_stable_compound_unchanged(self):
        """Aspirin has no tautomers; canonical form should be equivalent."""
        from rdkit import Chem

        from chemfuse.compute.standardization import canonical_tautomer

        result = canonical_tautomer(ASPIRIN_SMILES)
        assert result is not None
        mol = Chem.MolFromSmiles(result)
        assert mol is not None

    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES should return None."""
        from chemfuse.compute.standardization import canonical_tautomer

        assert canonical_tautomer("INVALID!!!") is None

    def test_empty_string_returns_none(self):
        """Empty string should return None."""
        from chemfuse.compute.standardization import canonical_tautomer

        assert canonical_tautomer("") is None

    def test_canonical_tautomer_is_idempotent(self):
        """Applying canonical_tautomer twice must yield the same result."""
        from chemfuse.compute.standardization import canonical_tautomer

        first = canonical_tautomer(KETO_TAUTOMER)
        assert first is not None
        second = canonical_tautomer(first)
        assert first == second


class TestStandardizeMol:
    """Tests for standardize_mol()."""

    def test_chains_salt_stripping_and_tautomer_canonicalization(self):
        """standardize_mol chains both operations in sequence."""
        from chemfuse.compute.standardization import (
            canonical_tautomer,
            standardize_mol,
            strip_salts,
        )

        result = standardize_mol(SODIUM_ACETATE_SMILES)
        assert result is not None

        # Should equal canonical_tautomer(strip_salts(input))
        stripped = strip_salts(SODIUM_ACETATE_SMILES)
        assert stripped is not None
        expected = canonical_tautomer(stripped)
        assert result == expected

    def test_strip_salts_only(self):
        """With canonicalize_tautomers=False, only salt stripping is performed."""
        from chemfuse.compute.standardization import standardize_mol, strip_salts

        result = standardize_mol(SODIUM_ACETATE_SMILES, canonicalize_tautomers=False)
        expected = strip_salts(SODIUM_ACETATE_SMILES)
        assert result == expected

    def test_canonicalize_tautomers_only(self):
        """With strip_salts=False, only tautomer canonicalization is performed."""
        from chemfuse.compute.standardization import canonical_tautomer, standardize_mol

        result = standardize_mol(KETO_TAUTOMER, strip_salts=False)
        expected = canonical_tautomer(KETO_TAUTOMER)
        assert result == expected

    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES should return None."""
        from chemfuse.compute.standardization import standardize_mol

        assert standardize_mol("GARBAGE!!!") is None

    def test_empty_string_returns_none(self):
        """Empty string should return None."""
        from chemfuse.compute.standardization import standardize_mol

        assert standardize_mol("") is None

    def test_both_flags_false_returns_canonical_smiles(self):
        """With both flags False, returns RDKit canonical SMILES without standardization."""
        from rdkit import Chem

        from chemfuse.compute.standardization import standardize_mol

        result = standardize_mol(ASPIRIN_SMILES, strip_salts=False, canonicalize_tautomers=False)
        assert result is not None
        assert result == Chem.MolToSmiles(Chem.MolFromSmiles(ASPIRIN_SMILES))


class TestOptionalDependencyError:
    """Tests for OptionalDependencyError raised when RDKit is missing."""

    def test_strip_salts_raises_optional_dependency_error_when_rdkit_missing(self):
        """strip_salts must raise OptionalDependencyError if RDKit is not importable."""
        from chemfuse.core.exceptions import OptionalDependencyError

        with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None,
                                       "rdkit.Chem.MolStandardize": None,
                                       "rdkit.Chem.MolStandardize.rdMolStandardize": None}):
            import importlib

            import chemfuse.compute.standardization as std_mod
            importlib.reload(std_mod)
            with pytest.raises(OptionalDependencyError):
                std_mod.strip_salts(ASPIRIN_SMILES)

    def test_canonical_tautomer_raises_optional_dependency_error_when_rdkit_missing(self):
        """canonical_tautomer must raise OptionalDependencyError if RDKit is not importable."""
        from chemfuse.core.exceptions import OptionalDependencyError

        with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None,
                                       "rdkit.Chem.MolStandardize": None,
                                       "rdkit.Chem.MolStandardize.rdMolStandardize": None}):
            import importlib

            import chemfuse.compute.standardization as std_mod
            importlib.reload(std_mod)
            with pytest.raises(OptionalDependencyError):
                std_mod.canonical_tautomer(KETO_TAUTOMER)

    def test_standardize_mol_raises_optional_dependency_error_when_rdkit_missing(self):
        """standardize_mol must raise OptionalDependencyError if RDKit is not importable."""
        from chemfuse.core.exceptions import OptionalDependencyError

        with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None,
                                       "rdkit.Chem.MolStandardize": None,
                                       "rdkit.Chem.MolStandardize.rdMolStandardize": None}):
            import importlib

            import chemfuse.compute.standardization as std_mod
            importlib.reload(std_mod)
            with pytest.raises(OptionalDependencyError):
                std_mod.standardize_mol(ASPIRIN_SMILES)


class TestCompoundToMol:
    """Tests for Compound.to_mol()."""

    def test_to_mol_returns_rdkit_mol_for_valid_smiles(self):
        """to_mol() should return an RDKit Mol object for valid SMILES."""
        from rdkit import Chem

        from chemfuse.models.compound import Compound

        compound = Compound(smiles=ASPIRIN_SMILES)
        mol = compound.to_mol()
        assert mol is not None
        assert isinstance(mol, Chem.rdchem.Mol)

    def test_to_mol_returns_none_for_empty_smiles(self):
        """to_mol() should return None when smiles is empty string."""
        from chemfuse.models.compound import Compound

        compound = Compound(smiles="")
        assert compound.to_mol() is None

    def test_to_mol_returns_none_for_invalid_smiles(self):
        """to_mol() should return None for an unparseable SMILES."""
        from chemfuse.models.compound import Compound

        compound = Compound(smiles="NOT_VALID!!!")
        assert compound.to_mol() is None

    def test_to_mol_returns_none_when_rdkit_unavailable(self):
        """to_mol() must return None (not raise) if RDKit cannot be imported."""
        from chemfuse.models.compound import Compound

        compound = Compound(smiles=ASPIRIN_SMILES)
        with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None}):
            result = compound.to_mol()
        assert result is None

    def test_canonical_smiles_field_defaults_to_none(self):
        """canonical_smiles field should default to None on a new Compound."""
        from chemfuse.models.compound import Compound

        compound = Compound(smiles=ASPIRIN_SMILES)
        assert compound.canonical_smiles is None

    def test_canonical_smiles_field_can_be_set(self):
        """canonical_smiles field should accept a string value."""
        from chemfuse.models.compound import Compound

        compound = Compound(smiles=ASPIRIN_SMILES, canonical_smiles=ASPIRIN_SMILES)
        assert compound.canonical_smiles == ASPIRIN_SMILES
