"""Compound standardization: salt stripping, tautomer canonicalization.

Provides structure normalization utilities to ensure consistent SMILES
representations for deduplication and descriptor computation.
"""
from __future__ import annotations


def strip_salts(smiles: str) -> str | None:
    """Return the largest organic fragment (strips counterions and salts).

    Uses RDKit's LargestFragmentChooser to select the heaviest fragment,
    removing inorganic counterions such as Na+, Cl-, K+, etc.

    Args:
        smiles: Input SMILES string (may contain dot-disconnected fragments).

    Returns:
        Canonical SMILES of the largest fragment, or None if input is
        invalid or empty.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:
        from chemfuse.core.exceptions import OptionalDependencyError
        raise OptionalDependencyError("rdkit", "compute") from exc

    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    chooser = rdMolStandardize.LargestFragmentChooser()
    largest = chooser.choose(mol)
    if largest is None:
        return None

    # Neutralize formal charges (e.g. acetate [O-] -> OH, ammonium [NH4+] -> NH3)
    uncharger = rdMolStandardize.Uncharger()
    largest = uncharger.uncharge(largest)
    if largest is None:
        return None

    return Chem.MolToSmiles(largest)


def canonical_tautomer(smiles: str) -> str | None:
    """Return the canonical tautomer SMILES.

    Uses RDKit's TautomerEnumerator to enumerate all tautomers and return
    the canonical (lowest-energy) representative, ensuring consistent
    representation across keto/enol, lactam/lactim, etc.

    Args:
        smiles: Input SMILES string.

    Returns:
        Canonical tautomer SMILES, or None if input is invalid or empty.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:
        from chemfuse.core.exceptions import OptionalDependencyError
        raise OptionalDependencyError("rdkit", "compute") from exc

    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    enumerator = rdMolStandardize.TautomerEnumerator()
    canonical = enumerator.Canonicalize(mol)
    if canonical is None:
        return None

    return Chem.MolToSmiles(canonical)


def standardize_mol(
    smiles: str,
    strip_salts: bool = True,
    canonicalize_tautomers: bool = True,
) -> str | None:
    """Standardize a SMILES string by stripping salts and/or canonicalizing tautomers.

    Chains both operations: salt stripping is applied first, then tautomer
    canonicalization on the resulting fragment.

    Args:
        smiles: Input SMILES string.
        strip_salts: If True, remove counterions and select the largest fragment.
        canonicalize_tautomers: If True, normalize to canonical tautomer form.

    Returns:
        Standardized canonical SMILES, or None if input is invalid or empty.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:
        from chemfuse.core.exceptions import OptionalDependencyError
        raise OptionalDependencyError("rdkit", "compute") from exc

    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    if strip_salts:
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
        if mol is None:
            return None
        # Neutralize formal charges after fragment selection
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        if mol is None:
            return None

    if canonicalize_tautomers:
        enumerator = rdMolStandardize.TautomerEnumerator()
        mol = enumerator.Canonicalize(mol)
        if mol is None:
            return None

    return Chem.MolToSmiles(mol)
