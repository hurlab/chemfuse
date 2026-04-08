"""Synthetic accessibility, fraction sp3, and natural product-likeness scores.

Provides three molecular property functions:
    - synthetic_accessibility: SA_Score via RDKit Contrib (1=easy, 10=hard)
    - fraction_sp3: Fraction of sp3 carbons (Fsp3), 0.0 to 1.0
    - np_likeness: NP-likeness score via RDKit Contrib (optional)

All functions return None for invalid SMILES or when RDKit is unavailable.
"""

from __future__ import annotations

import logging
import os
import sys
from typing import TYPE_CHECKING

from chemfuse.core.exceptions import OptionalDependencyError

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# RDKit availability check
# ---------------------------------------------------------------------------

try:
    from rdkit import Chem, RDConfig  # type: ignore[import-not-found]
    from rdkit.Chem import rdMolDescriptors  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

# ---------------------------------------------------------------------------
# SA_Score loader (lazy, module-level singleton)
# ---------------------------------------------------------------------------

_sascorer = None
_sascorer_loaded = False


def _load_sascorer() -> object | None:
    """Load sascorer from RDKit Contrib directory.

    Returns the sascorer module or None if unavailable.
    """
    global _sascorer, _sascorer_loaded  # noqa: PLW0603
    if _sascorer_loaded:
        return _sascorer

    _sascorer_loaded = True
    if not _RDKIT_AVAILABLE:
        return None

    try:
        sa_score_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_score_dir not in sys.path:
            sys.path.append(sa_score_dir)
        import sascorer  # type: ignore[import-not-found]

        _sascorer = sascorer
    except ImportError:
        logger.debug("RDKit SA_Score contrib not available; SA score will be None.")
        _sascorer = None

    return _sascorer


# ---------------------------------------------------------------------------
# NP_Score loader (lazy, module-level singleton)
# ---------------------------------------------------------------------------

_npscorer = None
_npscorer_loaded = False
_np_model = None


def _load_npscorer() -> object | None:
    """Load npscorer from RDKit Contrib directory.

    Returns the npscorer module or None if unavailable.
    """
    global _npscorer, _npscorer_loaded, _np_model  # noqa: PLW0603
    if _npscorer_loaded:
        return _npscorer

    _npscorer_loaded = True
    if not _RDKIT_AVAILABLE:
        return None

    try:
        np_score_dir = os.path.join(RDConfig.RDContribDir, "NP_Score")
        if np_score_dir not in sys.path:
            sys.path.append(np_score_dir)
        import npscorer  # type: ignore[import-not-found]

        _npscorer = npscorer
        _np_model = npscorer.readNPModel()
    except (ImportError, OSError) as exc:
        logger.debug("RDKit NP_Score contrib not available: %s", exc)
        _npscorer = None

    return _npscorer


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def synthetic_accessibility(smiles: str) -> float | None:
    """Compute the synthetic accessibility (SA) score for a molecule.

    Uses RDKit's SA_Score from the Contrib directory. The score ranges from
    1.0 (trivially easy to synthesize) to 10.0 (extremely difficult).

    Args:
        smiles: SMILES string of the molecule.

    Returns:
        SA score as a float in [1.0, 10.0], or None if the SMILES is invalid
        or if RDKit / SA_Score is not available.

    Raises:
        OptionalDependencyError: If RDKit is not installed at all.
    """
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", extra="chem")

    if not smiles or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    sascorer = _load_sascorer()
    if sascorer is None:
        return None

    try:
        score: float = sascorer.calculateScore(mol)
        return score
    except Exception as exc:
        logger.warning("SA score computation failed for %r: %s", smiles, exc)
        return None


def fraction_sp3(smiles: str) -> float | None:
    """Compute the fraction of sp3 carbons (Fsp3) for a molecule.

    Fsp3 = (number of sp3 carbons) / (total number of carbons).
    Values range from 0.0 (entirely aromatic/sp2) to 1.0 (fully saturated).

    Args:
        smiles: SMILES string of the molecule.

    Returns:
        Fsp3 as a float in [0.0, 1.0], or None if the SMILES is invalid
        or if RDKit is not available.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", extra="chem")

    if not smiles or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        fsp3: float = rdMolDescriptors.CalcFractionCSP3(mol)
        return fsp3
    except Exception as exc:
        logger.warning("Fsp3 computation failed for %r: %s", smiles, exc)
        return None


def np_likeness(smiles: str) -> float | None:
    """Compute the natural product (NP) likeness score for a molecule.

    Uses RDKit's NP_Score from the Contrib directory. Scores typically range
    from approximately -5 (not NP-like) to +5 (very NP-like).

    This function is fully optional: if the NP_Score contrib is unavailable,
    it returns None silently without raising an error.

    Args:
        smiles: SMILES string of the molecule.

    Returns:
        NP-likeness score as a float, or None if the SMILES is invalid or if
        the NP_Score module is not available.

    Raises:
        OptionalDependencyError: If RDKit is not installed at all.
    """
    if not _RDKIT_AVAILABLE:
        raise OptionalDependencyError("rdkit", extra="chem")

    if not smiles or not smiles.strip():
        return None

    npscorer = _load_npscorer()
    if npscorer is None or _np_model is None:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        score: float = npscorer.scoreMol(mol, _np_model)
        return score
    except Exception as exc:
        logger.warning("NP-likeness computation failed for %r: %s", smiles, exc)
        return None
