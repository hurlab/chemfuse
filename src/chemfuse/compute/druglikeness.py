"""Drug-likeness filters for ChemFuse.

Implements the five standard drug-likeness filters:
    - Lipinski Rule of Five (Lipinski et al., 1997)
    - Veber Rules (Veber et al., 2002)
    - Ghose Filter (Ghose et al., 1999)
    - Egan Filter (Egan et al., 2000)
    - Muegge Filter (Muegge et al., 2001)

Plus two additional analyses:
    - PAINS Filter (Baell & Holloway, 2010, J Med Chem 53:2719)
    - QED Score (Bickerton et al., 2012, Nature Chem 4:90)

Filters work WITHOUT RDKit by using pre-computed properties (e.g. from PubChem).
When RDKit IS installed and a SMILES string is provided, properties are computed
locally, filling in any None values from the fetched properties.

PAINS and QED require RDKit and a valid SMILES or molecule object.
"""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.models.prediction import DrugLikeness, FilterResult

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem  # type: ignore[import-not-found]
    from rdkit.Chem import (
        Descriptors,  # type: ignore[import-not-found]
        rdMolDescriptors,  # type: ignore[import-not-found]
    )
    from rdkit.Chem.FilterCatalog import (  # type: ignore[import-not-found]
        FilterCatalog,
        FilterCatalogParams,
    )
    from rdkit.Chem.QED import qed as _rdkit_qed  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True

    # Build the PAINS catalog once at module load; it is expensive to initialise.
    _pains_params = FilterCatalogParams()
    _pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    _PAINS_CATALOG: FilterCatalog = FilterCatalog(_pains_params)

except ImportError:
    _RDKIT_AVAILABLE = False
    _PAINS_CATALOG = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fill_from_rdkit(
    smiles: str,
    props: dict[str, Any],
) -> dict[str, Any]:
    """Fill None values in *props* using RDKit when SMILES is provided.

    Only called when ``_RDKIT_AVAILABLE`` is True.  Mutates *props* in place
    and returns it.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return props

    if props.get("mw") is None:
        props["mw"] = Descriptors.MolWt(mol)
    if props.get("logp") is None:
        props["logp"] = Descriptors.MolLogP(mol)
    if props.get("hbd") is None:
        props["hbd"] = rdMolDescriptors.CalcNumHBD(mol)
    if props.get("hba") is None:
        props["hba"] = rdMolDescriptors.CalcNumHBA(mol)
    if props.get("tpsa") is None:
        props["tpsa"] = Descriptors.TPSA(mol)
    if props.get("rotatable_bonds") is None:
        props["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if props.get("molar_refractivity") is None:
        props["molar_refractivity"] = Descriptors.MolMR(mol)
    if props.get("heavy_atom_count") is None:
        props["heavy_atom_count"] = mol.GetNumHeavyAtoms()
    if props.get("total_atoms") is None:
        # Ghose 1999 uses total atom count including implicit hydrogens.
        # AddHs() adds explicit H then GetNumAtoms() gives the full count.
        mol_with_h = Chem.AddHs(mol)
        props["total_atoms"] = mol_with_h.GetNumAtoms()
    if props.get("ring_count") is None:
        props["ring_count"] = rdMolDescriptors.CalcNumRings(mol)
    if props.get("heteroatom_count") is None:
        props["heteroatom_count"] = rdMolDescriptors.CalcNumHeteroatoms(mol)
    if props.get("carbon_count") is None:
        props["carbon_count"] = sum(
            1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6
        )

    return props


def _build_props(
    smiles: str | None,
    molecular_weight: float | None,
    logp: float | None,
    hbd: int | None,
    hba: int | None,
    tpsa: float | None,
    rotatable_bonds: int | None,
    molar_refractivity: float | None = None,
    total_atoms: int | None = None,
    heavy_atom_count: int | None = None,
    ring_count: int | None = None,
    heteroatom_count: int | None = None,
    carbon_count: int | None = None,
) -> dict[str, Any]:
    """Collect provided properties and optionally fill gaps via RDKit."""
    props: dict[str, Any] = {
        "mw": molecular_weight,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "tpsa": tpsa,
        "rotatable_bonds": rotatable_bonds,
        "molar_refractivity": molar_refractivity,
        # Keep total_atoms and heavy_atom_count separate.
        # total_atoms (all atoms including H) is filled by RDKit when SMILES is
        # provided. heavy_atom_count from PubChem is NOT a valid proxy because
        # Ghose's atom count includes hydrogens.
        "total_atoms": total_atoms,
        "heavy_atom_count": heavy_atom_count,
        "ring_count": ring_count,
        "heteroatom_count": heteroatom_count,
        "carbon_count": carbon_count,
    }

    if smiles and _RDKIT_AVAILABLE:
        props = _fill_from_rdkit(smiles, props)

    return props


# ---------------------------------------------------------------------------
# Individual filter functions
# ---------------------------------------------------------------------------

def lipinski_filter(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> FilterResult:
    """Apply Lipinski's Rule of Five.

    Criteria: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10.
    Passes with at most 1 violation (original Lipinski rule).

    Args:
        properties: Dict with keys 'molecular_weight', 'xlogp', 'hbd_count',
            'hba_count' (uses PubChem field names or generic aliases).
        smiles: Optional SMILES; used to compute missing properties via RDKit.

    Returns:
        FilterResult with pass_filter, violations, and details.
    """
    mw = _extract(properties, "molecular_weight", "mw")
    logp = _extract(properties, "xlogp", "logp")
    hbd = _extract_int(properties, "hbd_count", "hbd")
    hba = _extract_int(properties, "hba_count", "hba")

    props = _build_props(
        smiles, mw, logp, hbd, hba, None, None,
    )
    mw = props["mw"]
    logp = props["logp"]
    hbd = props["hbd"]
    hba = props["hba"]

    # Check for completely missing data
    if mw is None and logp is None and hbd is None and hba is None:
        return FilterResult(
            pass_filter=False,
            violations=["Missing required properties: molecular_weight, xlogp, hbd_count, hba_count"],
            details={},
        )

    violations: list[str] = []

    if mw is not None and mw > 500:
        violations.append(f"MW {mw:.2f} exceeds threshold 500")
    if logp is not None and logp > 5:
        violations.append(f"LogP {logp:.2f} exceeds threshold 5")
    if hbd is not None and hbd > 5:
        violations.append(f"HBD {hbd} exceeds threshold 5")
    if hba is not None and hba > 10:
        violations.append(f"HBA {hba} exceeds threshold 10")

    # Lipinski allows 1 violation
    pass_filter = len(violations) <= 1

    details: dict[str, object] = {
        "molecular_weight": {"value": mw, "threshold": "<=500"},
        "logp": {"value": logp, "threshold": "<=5"},
        "hbd": {"value": hbd, "threshold": "<=5"},
        "hba": {"value": hba, "threshold": "<=10"},
        "violations_count": len(violations),
        "allowed_violations": 1,
    }

    return FilterResult(pass_filter=pass_filter, violations=violations, details=details)


def veber_filter(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> FilterResult:
    """Apply Veber oral bioavailability rules.

    Criteria: TPSA <= 140, Rotatable Bonds <= 10.
    No violations allowed.

    Args:
        properties: Dict with keys 'tpsa' and 'rotatable_bonds' (or aliases).
        smiles: Optional SMILES; used to compute missing properties via RDKit.

    Returns:
        FilterResult with pass_filter, violations, and details.
    """
    tpsa = _extract(properties, "tpsa")
    rotatable_bonds = _extract_int(properties, "rotatable_bonds", "rotatable_bond_count")

    props = _build_props(smiles, None, None, None, None, tpsa, rotatable_bonds)
    tpsa = props["tpsa"]
    rotatable_bonds = props["rotatable_bonds"]

    if tpsa is None and rotatable_bonds is None:
        return FilterResult(
            pass_filter=False,
            violations=["Missing required properties: tpsa, rotatable_bonds"],
            details={},
        )

    violations: list[str] = []

    if tpsa is not None and tpsa > 140:
        violations.append(f"TPSA {tpsa:.2f} exceeds threshold 140")
    if rotatable_bonds is not None and rotatable_bonds > 10:
        violations.append(f"RotBonds {rotatable_bonds} exceeds threshold 10")

    details: dict[str, object] = {
        "tpsa": {"value": tpsa, "threshold": "<=140"},
        "rotatable_bonds": {"value": rotatable_bonds, "threshold": "<=10"},
    }

    return FilterResult(
        pass_filter=len(violations) == 0,
        violations=violations,
        details=details,
    )


def ghose_filter(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> FilterResult:
    """Apply Ghose drug-likeness filter (Ghose et al., 1999).

    Criteria:
        160 <= MW <= 480
        -0.4 <= LogP <= 5.6
        40 <= Molar Refractivity <= 130
        20 <= Heavy Atom Count (total atoms) <= 70

    No violations allowed.

    Args:
        properties: Dict with relevant property keys.
        smiles: Optional SMILES; used to compute missing properties via RDKit.

    Returns:
        FilterResult with pass_filter, violations, and details.
    """
    mw = _extract(properties, "molecular_weight", "mw")
    logp = _extract(properties, "xlogp", "logp")
    mr = _extract(properties, "molar_refractivity", "mr")
    # total_atoms in Ghose is ALL atoms including hydrogens.  PubChem's
    # heavy_atom_count only counts non-H atoms, so it cannot be used as a
    # proxy.  We leave total_atoms=None here; _fill_from_rdkit will set it to
    # the correct value (heavy + H) when SMILES is available.  When SMILES is
    # not available the atoms check is skipped (value stays None).
    heavy_atom_count = _extract_int(properties, "heavy_atom_count")

    props = _build_props(
        smiles, mw, logp, None, None, None, None,
        molar_refractivity=mr,
        total_atoms=None,
        heavy_atom_count=heavy_atom_count,
    )
    mw = props["mw"]
    logp = props["logp"]
    mr = props["molar_refractivity"]
    total_atoms = props["total_atoms"]

    if mw is None and logp is None and mr is None and total_atoms is None:
        return FilterResult(
            pass_filter=False,
            violations=["Insufficient molecular data"],
            details={},
        )

    violations: list[str] = []

    if mw is not None and (mw < 160 or mw > 480):
        violations.append(f"MW {mw:.2f} not in [160, 480]")
    if logp is not None and (logp < -0.4 or logp > 5.6):
        violations.append(f"LogP {logp:.2f} not in [-0.4, 5.6]")
    if mr is not None and (mr < 40 or mr > 130):
        violations.append(f"MR {mr:.2f} not in [40, 130]")
    if total_atoms is not None and (total_atoms < 20 or total_atoms > 70):
        violations.append(f"Atoms {total_atoms} not in [20, 70]")

    details: dict[str, object] = {
        "molecular_weight": {"value": mw, "threshold": "[160, 480]"},
        "logp": {"value": logp, "threshold": "[-0.4, 5.6]"},
        "molar_refractivity": {"value": mr, "threshold": "[40, 130]"},
        "total_atoms": {"value": total_atoms, "threshold": "[20, 70]"},
    }

    return FilterResult(
        pass_filter=len(violations) == 0,
        violations=violations,
        details=details,
    )


def egan_filter(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> FilterResult:
    """Apply Egan absorption filter (Egan et al., 2000).

    Criteria:
        -1.0 <= LogP <= 5.88
        TPSA <= 131.6

    No violations allowed.

    Args:
        properties: Dict with 'xlogp'/'logp' and 'tpsa' keys.
        smiles: Optional SMILES; used to compute missing properties via RDKit.

    Returns:
        FilterResult with pass_filter, violations, and details.
    """
    logp = _extract(properties, "xlogp", "logp")
    tpsa = _extract(properties, "tpsa")

    props = _build_props(smiles, None, logp, None, None, tpsa, None)
    logp = props["logp"]
    tpsa = props["tpsa"]

    if logp is None and tpsa is None:
        return FilterResult(
            pass_filter=False,
            violations=["Insufficient molecular data"],
            details={},
        )

    violations: list[str] = []

    if logp is not None and (logp < -1.0 or logp > 5.88):
        violations.append(f"LogP {logp:.2f} not in [-1.0, 5.88]")
    if tpsa is not None and tpsa > 131.6:
        violations.append(f"TPSA {tpsa:.2f} exceeds threshold 131.6")

    details: dict[str, object] = {
        "logp": {"value": logp, "threshold": "[-1.0, 5.88]"},
        "tpsa": {"value": tpsa, "threshold": "<=131.6"},
    }

    return FilterResult(
        pass_filter=len(violations) == 0,
        violations=violations,
        details=details,
    )


def muegge_filter(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> FilterResult:
    """Apply Muegge drug-likeness filter (Muegge et al., 2001).

    Criteria:
        160 <= MW <= 600
        -2 <= LogP <= 5
        TPSA <= 150
        Rings <= 7
        HBA <= 10
        HBD <= 5
        Rotatable Bonds <= 15
        Carbons > 4
        Heteroatoms > 1

    No violations allowed.

    Args:
        properties: Dict with relevant property keys.
        smiles: Optional SMILES; used to compute missing properties via RDKit.

    Returns:
        FilterResult with pass_filter, violations, and details.
    """
    mw = _extract(properties, "molecular_weight", "mw")
    logp = _extract(properties, "xlogp", "logp")
    tpsa = _extract(properties, "tpsa")
    hbd = _extract_int(properties, "hbd_count", "hbd")
    hba = _extract_int(properties, "hba_count", "hba")
    rotatable_bonds = _extract_int(properties, "rotatable_bonds", "rotatable_bond_count")
    ring_count = _extract_int(properties, "ring_count", "num_rings")
    heteroatom_count = _extract_int(properties, "heteroatom_count", "num_heteroatoms")
    carbon_count = _extract_int(properties, "carbon_count", "num_carbons")

    props = _build_props(
        smiles, mw, logp, hbd, hba, tpsa, rotatable_bonds,
        ring_count=ring_count,
        heteroatom_count=heteroatom_count,
        carbon_count=carbon_count,
    )
    mw = props["mw"]
    logp = props["logp"]
    tpsa = props["tpsa"]
    hbd = props["hbd"]
    hba = props["hba"]
    rotatable_bonds = props["rotatable_bonds"]
    ring_count = props["ring_count"]
    heteroatom_count = props["heteroatom_count"]
    carbon_count = props["carbon_count"]

    if all(
        v is None
        for v in [mw, logp, tpsa, hbd, hba, rotatable_bonds, ring_count, heteroatom_count, carbon_count]
    ):
        return FilterResult(
            pass_filter=False,
            violations=["Insufficient molecular data"],
            details={},
        )

    violations: list[str] = []

    if mw is not None and (mw < 160 or mw > 600):
        violations.append(f"MW {mw:.2f} not in [160, 600]")
    if logp is not None and (logp < -2 or logp > 5):
        violations.append(f"LogP {logp:.2f} not in [-2, 5]")
    if tpsa is not None and tpsa > 150:
        violations.append(f"TPSA {tpsa:.2f} exceeds threshold 150")
    if hbd is not None and hbd > 5:
        violations.append(f"HBD {hbd} exceeds threshold 5")
    if hba is not None and hba > 10:
        violations.append(f"HBA {hba} exceeds threshold 10")
    if ring_count is not None and ring_count > 7:
        violations.append(f"Rings {ring_count} exceeds threshold 7")
    if rotatable_bonds is not None and rotatable_bonds > 15:
        violations.append(f"RotBonds {rotatable_bonds} exceeds threshold 15")
    if heteroatom_count is not None and heteroatom_count <= 1:
        violations.append(f"Heteroatoms {heteroatom_count} not > 1")
    if carbon_count is not None and carbon_count <= 4:
        violations.append(f"Carbons {carbon_count} not > 4")

    details: dict[str, object] = {
        "molecular_weight": {"value": mw, "threshold": "[160, 600]"},
        "logp": {"value": logp, "threshold": "[-2, 5]"},
        "tpsa": {"value": tpsa, "threshold": "<=150"},
        "hbd": {"value": hbd, "threshold": "<=5"},
        "hba": {"value": hba, "threshold": "<=10"},
        "ring_count": {"value": ring_count, "threshold": "<=7"},
        "rotatable_bonds": {"value": rotatable_bonds, "threshold": "<=15"},
        "heteroatom_count": {"value": heteroatom_count, "threshold": ">1"},
        "carbon_count": {"value": carbon_count, "threshold": ">4"},
    }

    return FilterResult(
        pass_filter=len(violations) == 0,
        violations=violations,
        details=details,
    )


def pains_filter(
    smiles: str,
) -> FilterResult:
    """Apply the PAINS (Pan Assay Interference Compounds) substructure filter.

    Identifies compounds with substructures known to produce false positives in
    high-throughput screens (Baell & Holloway, 2010, J Med Chem 53:2719).

    Requires RDKit.  Raises ``ImportError`` when RDKit is not installed.

    Args:
        smiles: SMILES string for the compound to evaluate.

    Returns:
        FilterResult with:
            pass_filter=True  when no PAINS alert is matched,
            pass_filter=False when a PAINS substructure is found.
            violations: empty list or single-item list with the matched pattern
                name; violations count reflects 0 or 1 respectively.
            details: includes 'pains_match' key with the matched pattern name
                (or None when the compound passes).

    Raises:
        ImportError: when RDKit is not available.
        ValueError: when the SMILES string cannot be parsed.
    """
    if not _RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for pains_filter. "
            "Install it with: pip install rdkit"
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")

    entry = _PAINS_CATALOG.GetFirstMatch(mol)

    if entry is None:
        return FilterResult(
            pass_filter=True,
            violations=[],
            details={"pains_match": None},
        )

    pattern_name: str = entry.GetDescription()
    return FilterResult(
        pass_filter=False,
        violations=[f"PAINS alert matched: {pattern_name}"],
        details={"pains_match": pattern_name},
    )


# QED classification thresholds (Bickerton et al., 2012)
_QED_HIGH_THRESHOLD = 0.67
_QED_LOW_THRESHOLD = 0.33


def qed_score(
    smiles: str,
) -> dict[str, object]:
    """Compute the Quantitative Estimate of Drug-likeness (QED) score.

    Uses a desirability function over 8 molecular properties to produce a
    continuous score in [0, 1] (Bickerton et al., 2012, Nature Chem 4:90).

    Requires RDKit.  Raises ``ImportError`` when RDKit is not installed.

    Args:
        smiles: SMILES string for the compound to evaluate.

    Returns:
        Dict with:
            'qed' (float): score in [0, 1]; higher means more drug-like.
            'classification' (str): 'high' (>=0.67), 'medium' (>=0.33), or
                'low' (<0.33).

    Raises:
        ImportError: when RDKit is not available.
        ValueError: when the SMILES string cannot be parsed.
    """
    if not _RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for qed_score. "
            "Install it with: pip install rdkit"
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")

    score: float = _rdkit_qed(mol)

    if score >= _QED_HIGH_THRESHOLD:
        classification = "high"
    elif score >= _QED_LOW_THRESHOLD:
        classification = "medium"
    else:
        classification = "low"

    return {"qed": score, "classification": classification}


def check_drug_likeness(
    properties: dict[str, Any],
    smiles: str | None = None,
) -> DrugLikeness:
    """Run all drug-likeness filters and return a combined DrugLikeness result.

    Always runs the five core filters (Lipinski, Veber, Ghose, Egan, Muegge).
    When *smiles* is provided and RDKit is available, also runs the PAINS
    substructure filter and computes the QED score.

    Filters work with pre-computed properties (e.g. from PubChem) or with
    locally computed RDKit values when a SMILES string is provided.

    Args:
        properties: Dict with property values. Supported keys (PubChem names or
            generic aliases):
            - 'molecular_weight' / 'mw'
            - 'xlogp' / 'logp'
            - 'tpsa'
            - 'hbd_count' / 'hbd'
            - 'hba_count' / 'hba'
            - 'rotatable_bonds' / 'rotatable_bond_count'
            - 'heavy_atom_count' / 'total_atoms'
            - 'molar_refractivity' / 'mr'
            - 'ring_count' / 'num_rings'
            - 'heteroatom_count' / 'num_heteroatoms'
            - 'carbon_count' / 'num_carbons'
        smiles: Optional SMILES string; when provided with RDKit, fills missing
            properties automatically and enables PAINS + QED analysis.

    Returns:
        DrugLikeness object containing FilterResult for each of the five core
        filters, plus optional pains (FilterResult) and qed (dict) fields when
        SMILES and RDKit are available.
    """
    pains: FilterResult | None = None
    qed: dict[str, object] | None = None

    if smiles and _RDKIT_AVAILABLE:
        try:
            pains = pains_filter(smiles)
        except ValueError:
            logger.warning("PAINS filter skipped: could not parse SMILES %r", smiles)
        try:
            qed = qed_score(smiles)
        except ValueError:
            logger.warning("QED score skipped: could not parse SMILES %r", smiles)

    return DrugLikeness(
        lipinski=lipinski_filter(properties, smiles=smiles),
        veber=veber_filter(properties, smiles=smiles),
        ghose=ghose_filter(properties, smiles=smiles),
        egan=egan_filter(properties, smiles=smiles),
        muegge=muegge_filter(properties, smiles=smiles),
        pains=pains,
        qed=qed,
    )


# ---------------------------------------------------------------------------
# Utility extractors
# ---------------------------------------------------------------------------

def _extract(d: dict[str, Any], *keys: str) -> float | None:
    """Extract a float value from dict, trying multiple keys."""
    for key in keys:
        val = d.get(key)
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                pass
    return None


def _extract_int(d: dict[str, Any], *keys: str) -> int | None:
    """Extract an int value from dict, trying multiple keys."""
    for key in keys:
        val = d.get(key)
        if val is not None:
            try:
                return int(val)
            except (ValueError, TypeError):
                pass
    return None
