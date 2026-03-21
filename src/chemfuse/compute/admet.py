"""ADMET property prediction for ChemFuse.

Uses admet-ai (ML, Chemprop-based) when installed; falls back to rule-based
heuristics derived from physicochemical descriptors otherwise.
"""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.core.exceptions import ChemFuseError
from chemfuse.models.prediction import ADMETPrediction, ADMETProfile

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional dependency: admet-ai
# ---------------------------------------------------------------------------
admet_ai = None
try:
    import admet_ai  # type: ignore[import-not-found]

    _ADMET_AI_AVAILABLE = True
except ImportError:
    _ADMET_AI_AVAILABLE = False

# Optional dependency: RDKit (for rule-based fallback)
Chem = None
Descriptors = None
rdMolDescriptors = None
try:
    from rdkit import Chem  # type: ignore[import-not-found]
    from rdkit.Chem import Descriptors, rdMolDescriptors  # type: ignore[import-not-found]

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_smiles(smiles: str) -> Any:
    """Parse SMILES, raising ChemFuseError on failure."""
    if not smiles or not smiles.strip():
        raise ChemFuseError("SMILES string is empty.")
    if not _RDKIT_AVAILABLE:
        # Cannot validate without RDKit; return None and proceed with heuristics
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemFuseError(f"Invalid SMILES string: {smiles!r}")
    return mol


def _make_prediction(
    property_name: str,
    value: float | str,
    unit: str | None = None,
    confidence: float | None = None,
    method: str = "rule-based",
    category: str | None = None,
) -> ADMETPrediction:
    return ADMETPrediction(
        property_name=property_name,
        value=value,
        unit=unit,
        confidence=confidence,
        method=method,
        category=category,
    )


# ---------------------------------------------------------------------------
# Rule-based fallback predictions
# ---------------------------------------------------------------------------

def _rule_solubility(mol: Any) -> ADMETPrediction:
    """ESOL solubility: log(S) = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RotBonds - 0.74*AP."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    # Aromatic proportion
    aromatic_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = aromatic_atoms / heavy_atoms if heavy_atoms > 0 else 0.0

    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rot_bonds - 0.74 * ap

    if log_s >= -2:
        category = "high"
    elif log_s >= -4:
        category = "medium"
    else:
        category = "low"

    return _make_prediction(
        "solubility", round(log_s, 3), unit="log mol/L",
        confidence=0.7, method="rule-based", category=category,
    )


def _rule_gi_absorption(mol: Any) -> ADMETPrediction:
    """GI absorption: high if TPSA < 140 and LogP > -1."""
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    absorbed = tpsa < 140 and logp > -1
    category = "high" if (absorbed and tpsa < 80 and 0 < logp < 4) else ("medium" if absorbed else "low")
    return _make_prediction(
        "gi_absorption", 1.0 if absorbed else 0.0,
        confidence=0.65, method="rule-based", category=category,
    )


def _rule_bbb(mol: Any) -> ADMETPrediction:
    """BBB permeability: yes if TPSA < 90 and MW < 400."""
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    logp = Descriptors.MolLogP(mol)
    score = sum([mw < 400, tpsa < 90, hbd < 3, 1 <= logp <= 3]) / 4.0
    category = "high" if score >= 0.75 else ("medium" if score >= 0.5 else "low")
    return _make_prediction(
        "bbb_permeability", round(score, 3),
        confidence=0.6, method="rule-based", category=category,
    )


def _rule_cyp_inhibition(mol: Any, cyp: str) -> ADMETPrediction:
    """CYP inhibition risk based on aromatic ring count and nitrogen content."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    n_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    risk = sum([logp > 3, mw > 300, aromatic_rings >= 2, n_atoms >= 1]) / 4.0
    category = "high" if risk >= 0.75 else ("medium" if risk >= 0.5 else "low")
    return _make_prediction(
        cyp, round(risk, 3),
        confidence=0.55, method="rule-based", category=category,
    )


def _rule_pgp_substrate(mol: Any) -> ADMETPrediction:
    """P-gp substrate: based on MW > 400 or HBD > 3."""
    mw = Descriptors.MolWt(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    tpsa = Descriptors.TPSA(mol)
    score = sum([mw > 400, hbd > 3, tpsa > 100]) / 3.0
    category = "high" if score >= 0.67 else ("medium" if score >= 0.33 else "low")
    return _make_prediction(
        "pgp_substrate", round(score, 3),
        confidence=0.6, method="rule-based", category=category,
    )


def _rule_hepatotoxicity(mol: Any) -> ADMETPrediction:
    """Hepatotoxicity risk based on reactive group detection (logP, MW, TPSA)."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    risk = 0.0
    if logp > 3:
        risk += 0.3
    if mw > 500:
        risk += 0.2
    if tpsa < 25:
        risk += 0.2
    risk = min(risk, 1.0)
    category = "high" if risk >= 0.5 else ("medium" if risk >= 0.3 else "low")
    return _make_prediction(
        "hepatotoxicity", round(risk, 3),
        confidence=0.55, method="rule-based", category=category,
    )


def _rule_renal_clearance(mol: Any) -> ADMETPrediction:
    """Renal clearance based on MW and TPSA."""
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    if mw < 300 and tpsa > 60 and logp < 1:
        score, category = 0.8, "high"
    elif mw < 500 and logp < 3:
        score, category = 0.5, "medium"
    else:
        score, category = 0.2, "low"
    return _make_prediction(
        "renal_clearance", score,
        confidence=0.5, method="rule-based", category=category,
    )


def _rule_half_life(mol: Any) -> ADMETPrediction:
    """T1/2 heuristic: lipophilic compounds tend to have longer half-lives."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    # Simple empirical estimate in hours
    t12 = max(0.5, logp * 2.0 + mw / 200.0)
    t12 = round(min(t12, 48.0), 2)
    category = "high" if t12 > 12 else ("medium" if t12 > 4 else "low")
    return _make_prediction(
        "half_life", t12, unit="h",
        confidence=0.4, method="rule-based", category=category,
    )


def _rule_clearance(mol: Any) -> ADMETPrediction:
    """Clearance heuristic based on MW."""
    mw = Descriptors.MolWt(mol)
    # mL/min/kg rough estimate
    cl = round(max(1.0, 50.0 - mw / 20.0), 2)
    category = "high" if cl > 30 else ("medium" if cl > 10 else "low")
    return _make_prediction(
        "clearance", cl, unit="mL/min/kg",
        confidence=0.4, method="rule-based", category=category,
    )


def _rule_herg(mol: Any) -> ADMETPrediction:
    """hERG liability based on logP and molecular weight."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    n_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    # High logP + basic nitrogen = hERG risk
    risk = sum([logp > 3.7, mw > 300, n_atoms >= 1]) / 3.0
    category = "high" if risk >= 0.67 else ("medium" if risk >= 0.33 else "low")
    return _make_prediction(
        "herg_liability", round(risk, 3),
        confidence=0.55, method="rule-based", category=category,
    )


def _rule_ames(mol: Any) -> ADMETPrediction:
    """AMES mutagenicity: aromatic amines and nitro groups as risk factors."""
    smiles = Chem.MolToSmiles(mol)
    # Very simplified heuristic
    nitro = "N(=O)=O" in smiles or "[N+](=O)[O-]" in smiles
    aromatic_amine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
            aromatic_amine = True
            break
    score = (0.5 if nitro else 0.0) + (0.4 if aromatic_amine else 0.0)
    score = min(score, 1.0)
    category = "high" if score >= 0.5 else ("medium" if score >= 0.3 else "low")
    return _make_prediction(
        "ames_mutagenicity", round(score, 3),
        confidence=0.5, method="rule-based", category=category,
    )


def _rule_dili(mol: Any) -> ADMETPrediction:
    """DILI (Drug-Induced Liver Injury) risk."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    risk = sum([logp > 3, mw > 300]) / 2.0
    category = "high" if risk >= 0.75 else ("medium" if risk >= 0.5 else "low")
    return _make_prediction(
        "dili", round(risk, 3),
        confidence=0.5, method="rule-based", category=category,
    )


def _rule_lipophilicity(mol: Any) -> ADMETPrediction:
    """Log D (lipophilicity) approximation via MolLogP."""
    logp = Descriptors.MolLogP(mol)
    category = "high" if logp > 5 else ("medium" if logp > 1 else "low")
    return _make_prediction(
        "lipophilicity", round(logp, 3), unit="log D",
        confidence=0.75, method="rule-based", category=category,
    )


def _rule_caco2(mol: Any) -> ADMETPrediction:
    """Caco-2 permeability heuristic based on MW and TPSA."""
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    # nm/s estimate
    perm = max(0.0, 50.0 - tpsa / 3.0 - mw / 30.0)
    perm = round(min(perm, 50.0), 2)
    category = "high" if perm > 20 else ("medium" if perm > 5 else "low")
    return _make_prediction(
        "caco2_permeability", perm, unit="nm/s",
        confidence=0.55, method="rule-based", category=category,
    )


def _rule_hia(mol: Any) -> ADMETPrediction:
    """Human Intestinal Absorption."""
    tpsa = Descriptors.TPSA(mol)
    mw = Descriptors.MolWt(mol)
    # % absorbed estimate
    pct = max(0.0, min(100.0, 100.0 - tpsa / 2.0 - mw / 50.0))
    category = "high" if pct > 70 else ("medium" if pct > 30 else "low")
    return _make_prediction(
        "hia", round(pct, 1), unit="%",
        confidence=0.6, method="rule-based", category=category,
    )


def _rule_based_profile(smiles: str, mol: Any | None = None) -> dict[str, ADMETPrediction]:
    """Compute all rule-based ADMET predictions for a molecule."""
    if mol is None:
        if not _RDKIT_AVAILABLE:
            # Return placeholder predictions when RDKit is absent
            props = [
                "caco2_permeability", "hia", "bbb_permeability",
                "cyp1a2_inhibition", "cyp2c9_inhibition", "cyp2c19_inhibition",
                "cyp2d6_inhibition", "cyp3a4_inhibition", "half_life", "clearance",
                "herg_liability", "ames_mutagenicity", "dili",
                "solubility", "lipophilicity",
            ]
            return {
                p: _make_prediction(p, 0.0, confidence=None, method="rule-based")
                for p in props
            }
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemFuseError(f"Invalid SMILES string: {smiles!r}")

    preds: dict[str, ADMETPrediction] = {}
    preds["caco2_permeability"] = _rule_caco2(mol)
    preds["hia"] = _rule_hia(mol)
    preds["bbb_permeability"] = _rule_bbb(mol)
    preds["cyp1a2_inhibition"] = _rule_cyp_inhibition(mol, "cyp1a2_inhibition")
    preds["cyp2c9_inhibition"] = _rule_cyp_inhibition(mol, "cyp2c9_inhibition")
    preds["cyp2c19_inhibition"] = _rule_cyp_inhibition(mol, "cyp2c19_inhibition")
    preds["cyp2d6_inhibition"] = _rule_cyp_inhibition(mol, "cyp2d6_inhibition")
    preds["cyp3a4_inhibition"] = _rule_cyp_inhibition(mol, "cyp3a4_inhibition")
    preds["half_life"] = _rule_half_life(mol)
    preds["clearance"] = _rule_clearance(mol)
    preds["herg_liability"] = _rule_herg(mol)
    preds["ames_mutagenicity"] = _rule_ames(mol)
    preds["dili"] = _rule_dili(mol)
    preds["solubility"] = _rule_solubility(mol)
    preds["lipophilicity"] = _rule_lipophilicity(mol)
    return preds


def _overall_score(predictions: dict[str, ADMETPrediction]) -> float:
    """Compute an overall ADMET score in [0, 1] from category labels."""
    if not predictions:
        return 0.0
    cat_scores = {"low": 1.0, "medium": 0.5, "high": 0.0}
    # For solubility, hia, gi_absorption: higher is better
    positive_props = {"solubility", "gi_absorption", "hia", "caco2_permeability"}
    scores = []
    for name, pred in predictions.items():
        cat = pred.category
        if cat is None:
            continue
        if name in positive_props:
            score = {"high": 1.0, "medium": 0.5, "low": 0.0}.get(cat, 0.5)
        else:
            score = cat_scores.get(cat, 0.5)
        scores.append(score)
    return round(sum(scores) / len(scores), 3) if scores else 0.0


# ---------------------------------------------------------------------------
# ML-based predictions via admet-ai
# ---------------------------------------------------------------------------

def _ml_based_profile(smiles: str) -> dict[str, ADMETPrediction]:
    """Compute ADMET predictions using admet-ai ML models."""
    try:
        model = admet_ai.ADMETModel()  # type: ignore[name-defined]
        results: dict[str, Any] = model.predict(smiles=smiles)
    except Exception as exc:
        logger.warning("admet-ai prediction failed for %s: %s. Falling back to rule-based.", smiles, exc)
        return _rule_based_profile(smiles)

    preds: dict[str, ADMETPrediction] = {}
    for prop_name, value in results.items():
        if isinstance(value, dict):
            val = value.get("value", 0.0)
            conf = value.get("confidence", None)
        else:
            val = float(value) if value is not None else 0.0
            conf = None
        preds[str(prop_name)] = _make_prediction(
            str(prop_name), val, confidence=conf, method="ml"
        )
    return preds


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def predict_admet(smiles: str) -> ADMETProfile:
    """Predict ADMET properties for a single SMILES string.

    Uses admet-ai ML models when installed; falls back to rule-based
    heuristics otherwise.

    Args:
        smiles: Canonical SMILES string.

    Returns:
        ADMETProfile with all predictions populated.

    Raises:
        ChemFuseError: If the SMILES string is invalid.
    """
    mol = _validate_smiles(smiles)

    if _ADMET_AI_AVAILABLE:
        predictions = _ml_based_profile(smiles)
    else:
        predictions = _rule_based_profile(smiles, mol)

    score = _overall_score(predictions)
    return ADMETProfile(smiles=smiles, predictions=predictions, overall_score=score)


def predict_admet_batch(smiles_list: list[str]) -> list[ADMETProfile]:
    """Predict ADMET properties for a list of SMILES strings.

    Args:
        smiles_list: List of canonical SMILES strings.

    Returns:
        List of ADMETProfile objects, one per input SMILES.
    """
    results: list[ADMETProfile] = []
    for smiles in smiles_list:
        try:
            results.append(predict_admet(smiles))
        except ChemFuseError as exc:
            logger.warning("ADMET prediction failed for %r: %s", smiles, exc)
            results.append(ADMETProfile(smiles=smiles))
    return results
