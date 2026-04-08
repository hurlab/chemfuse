"""Natural language query interface for ChemFuse.

Parses natural language chemistry queries into structured API calls.
"""

from __future__ import annotations

import logging
import re
from typing import Any

logger = logging.getLogger(__name__)


def _is_smiles(token: str) -> bool:
    """Heuristically detect whether *token* looks like a SMILES string.

    A token is classified as SMILES if it:
    - Contains SMILES-only characters (=, #, [, ], @, +, /, \\)
    - Contains parentheses with adjacent alphanumeric characters (branch notation)
    - Contains lowercase aromatic atom symbols (c, n, o, p, s)
    - Contains digits adjacent to uppercase atom letters (ring closures)
    - Consists entirely of valid SMILES atom letters (C, N, O, P, S, B, F, I, Cl, Br)
      in a sequence that looks chemical (2+ chars, all uppercase SMILES atoms)

    Note: this intentionally works on the *original* (non-lowercased) token.
    Normalisation of SMILES to lowercase would break detection.
    """
    if not token:
        return False
    # Explicit SMILES-only characters
    if re.search(r"[=#\[\]@+/\\]", token):
        return True
    # Parentheses combined with letters or digits (branch notation)
    if re.search(r"\([A-Za-z0-9]|\)[A-Za-z0-9]|[A-Za-z0-9]\(", token):
        return True
    # Lowercase organic-subset atoms (aromatic atoms in SMILES): c, n, o, p, s
    # Match as standalone characters or adjacent to digits/capitals
    if re.search(r"(?<![A-Za-z])[cnops](?![A-Za-z])", token):
        return True
    # Digits embedded in a string that looks chemical (e.g. C1CCCCC1, N1, O2)
    if re.search(r"[A-Z]\d[A-Z]|[A-Z]\d$|\d[A-Z]", token):
        return True
    # Simple aliphatic SMILES: 2+ uppercase letters, each a valid SMILES atom symbol.
    # Covers CCO (ethanol), CC (ethane), CCCCC (pentane), etc.
    # Valid single-letter SMILES atoms: B, C, N, O, P, S, F, I
    if re.match(r"^[BCNOPSFIH][BCNOPSFIH]+$", token):
        return True
    return False


def _extract_smiles_tokens(text: str) -> list[str]:
    """Return all whitespace-separated tokens that look like SMILES."""
    return [tok for tok in text.split() if _is_smiles(tok)]


# ---------------------------------------------------------------------------
# Pattern definitions — all use re.IGNORECASE so they match normalised input.
# Each pattern captures the chemical entity as group(1) [and group(2) for compare].
# We apply each pattern against the ORIGINAL query (not lowercased) so that
# SMILES case is preserved in the captured groups.
# ---------------------------------------------------------------------------

_SIMILAR_KEYWORDS = re.compile(
    r"(?:similar\s+to|compounds?\s+like|find\s+similar(?:\s+to)?)\s+(.+)",
    re.IGNORECASE,
)

_ADMET_KEYWORDS = re.compile(
    r"(?:admet\s+(?:for|of|prediction\s+for)|predict\s+admet(?:\s+for)?)\s+(.+)",
    re.IGNORECASE,
)

# Handles both "is X drug-like" and "drug-likeness of X" and "check drug-likeness for X"
# Pattern A: keyword BEFORE entity  → drug-likeness of X / check drug-likeness for X
_DRUG_LIKE_PREFIX = re.compile(
    r"(?:check\s+drug[- ]like(?:ness)?(?:\s+for)?|drug[- ]like(?:ness)?\s+(?:for|of))\s+(.+)",
    re.IGNORECASE,
)
# Pattern B: entity BEFORE keyword  → is X drug-like
_DRUG_LIKE_SUFFIX = re.compile(
    r"(?:is\s+)(.+?)\s+drug[- ]like(?:ness)?",
    re.IGNORECASE,
)

_SEARCH_KEYWORDS = re.compile(
    r"(?:search\s+(?:for\s+)?|find\s+(?:compound\s+)?)(.+)",
    re.IGNORECASE,
)

_GET_KEYWORDS = re.compile(
    r"(?:get|look\s+up|retrieve|fetch)\s+(.+)",
    re.IGNORECASE,
)

_XREF_KEYWORDS = re.compile(
    r"(?:cross[- ]?reference|map|xref)\s+(.+)",
    re.IGNORECASE,
)

_COMPARE_KEYWORDS = re.compile(
    r"compare\s+(.+?)\s+and\s+(.+)",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Routing helpers
# ---------------------------------------------------------------------------

def _strip_trailing(query: str) -> str:
    """Collapse whitespace and strip trailing punctuation."""
    q = re.sub(r"\s+", " ", query.strip())
    return q.rstrip("?!")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# @MX:ANCHOR: ask() is the sole entry point for natural language queries
# @MX:REASON: Exported in __init__.py __all__; referenced externally
def ask(query: str) -> Any:
    """Execute a natural language chemistry query.

    Parses the query to determine the intended operation, then
    executes the corresponding ChemFuse API call.

    Supported query patterns:
    - "search for [compound]" / "find [compound]"
    - "get [identifier]" / "look up [identifier]"
    - "similar to [SMILES]" / "find compounds like [SMILES]"
    - "ADMET for [SMILES]" / "predict ADMET [SMILES]"
    - "is [SMILES] drug-like?" / "drug-likeness of [SMILES]"
    - "cross-reference [identifier]" / "map [identifier]"
    - "compare [SMILES1] and [SMILES2]"

    Args:
        query: Natural language query string.

    Returns:
        Query result (CompoundCollection, dict, etc.) depending on the query type.

    Raises:
        ValueError: If the query cannot be parsed.
    """
    if not query or not query.strip():
        raise ValueError("Query must not be empty.")

    # Preserve original casing for SMILES extraction; strip punctuation only.
    orig = _strip_trailing(query)

    # --- compare ---
    m = _COMPARE_KEYWORDS.search(orig)
    if m:
        return _handle_compare(m.group(1).strip(), m.group(2).strip(), orig)

    # --- find similar ---
    m = _SIMILAR_KEYWORDS.search(orig)
    if m:
        return _handle_similar(m.group(1).strip(), orig)

    # --- ADMET prediction ---
    m = _ADMET_KEYWORDS.search(orig)
    if m:
        return _handle_admet(m.group(1).strip(), orig)

    # --- drug-likeness (prefix: "drug-likeness of X") ---
    m = _DRUG_LIKE_PREFIX.search(orig)
    if m:
        return _handle_druglikeness(m.group(1).strip(), orig)

    # --- drug-likeness (suffix: "is X drug-like") ---
    m = _DRUG_LIKE_SUFFIX.search(orig)
    if m:
        return _handle_druglikeness(m.group(1).strip(), orig)

    # --- cross-reference / map ---
    m = _XREF_KEYWORDS.search(orig)
    if m:
        return _handle_xref(m.group(1).strip(), orig)

    # --- get / look up ---
    m = _GET_KEYWORDS.search(orig)
    if m:
        return _handle_get(m.group(1).strip(), orig)

    # --- search / find (must come after similar/get to avoid over-matching) ---
    m = _SEARCH_KEYWORDS.search(orig)
    if m:
        return _handle_search(m.group(1).strip(), orig)

    raise ValueError(
        f"Cannot parse query: {query!r}. "
        "Try 'search for aspirin', 'ADMET for CC(=O)O', or 'find similar to CC(=O)O'."
    )


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

def _handle_search(entity: str, original_query: str) -> Any:
    """Route a 'search for X' query to chemfuse.search()."""
    import chemfuse

    entity = entity.strip()
    if not entity:
        raise ValueError(f"No search term found in query: {original_query!r}")

    logger.debug("nlq.ask: search for %r", entity)

    if _is_smiles(entity):
        return chemfuse.search(entity, query_type="smiles")
    else:
        return chemfuse.search(entity, query_type="name")


def _handle_get(entity: str, original_query: str) -> Any:
    """Route a 'get X' query to chemfuse.get() or chemfuse.search()."""
    import chemfuse

    entity = entity.strip()
    if not entity:
        raise ValueError(f"No identifier found in query: {original_query!r}")

    logger.debug("nlq.ask: get %r", entity)

    # If purely numeric, treat as a PubChem CID
    if re.match(r"^\d+$", entity):
        return chemfuse.get(entity)

    if _is_smiles(entity):
        return chemfuse.search(entity, query_type="smiles")

    # Try name search
    return chemfuse.search(entity, query_type="name")


def _handle_similar(entity: str, original_query: str) -> Any:
    """Route a 'similar to X' query to chemfuse.find_similar()."""
    import chemfuse

    entity = entity.strip()
    smiles_tokens = _extract_smiles_tokens(entity)
    if smiles_tokens:
        smiles = smiles_tokens[0]
    elif _is_smiles(entity):
        smiles = entity
    else:
        # entity may be a compound name; look it up first
        coll = chemfuse.search(entity, query_type="name")
        if not coll:
            raise ValueError(
                f"Could not find SMILES for {entity!r} to run similarity search."
            )
        smiles = coll[0].smiles
        if not smiles:
            raise ValueError(
                f"Compound {entity!r} has no SMILES; cannot run similarity search."
            )

    logger.debug("nlq.ask: find_similar %r", smiles)
    return chemfuse.find_similar(smiles)


def _handle_admet(entity: str, original_query: str) -> Any:
    """Route an 'ADMET for X' query to chemfuse.compute.admet.predict_admet()."""
    from chemfuse.compute.admet import predict_admet

    entity = entity.strip()
    smiles_tokens = _extract_smiles_tokens(entity)

    if smiles_tokens:
        smiles = smiles_tokens[0]
    elif _is_smiles(entity):
        smiles = entity
    else:
        # Compound name: resolve to SMILES first
        import chemfuse
        coll = chemfuse.search(entity, query_type="name")
        if not coll:
            raise ValueError(f"Could not find SMILES for {entity!r} to predict ADMET.")
        smiles = coll[0].smiles
        if not smiles:
            raise ValueError(f"Compound {entity!r} has no SMILES; cannot predict ADMET.")

    logger.debug("nlq.ask: predict_admet %r", smiles)
    return predict_admet(smiles)


def _handle_druglikeness(entity: str, original_query: str) -> Any:
    """Route a 'drug-like X' query to compute.druglikeness.check_drug_likeness()."""
    from chemfuse.compute.druglikeness import check_drug_likeness

    entity = entity.strip()
    smiles_tokens = _extract_smiles_tokens(entity)

    if smiles_tokens:
        smiles = smiles_tokens[0]
    elif _is_smiles(entity):
        smiles = entity
    else:
        import chemfuse
        coll = chemfuse.search(entity, query_type="name")
        if not coll:
            raise ValueError(f"Could not find SMILES for {entity!r} to check drug-likeness.")
        smiles = coll[0].smiles
        if not smiles:
            raise ValueError(f"Compound {entity!r} has no SMILES; cannot check drug-likeness.")

    logger.debug("nlq.ask: check_drug_likeness %r", smiles)
    # check_drug_likeness accepts a properties dict; build one from SMILES when possible
    try:
        from chemfuse.compute.descriptors import compute_descriptors
        props = compute_descriptors(smiles)
    except Exception:  # noqa: BLE001
        props = {}
    return check_drug_likeness(props, smiles=smiles)


def _handle_xref(entity: str, original_query: str) -> Any:
    """Route a 'cross-reference X' query to chemfuse.map_identifiers()."""
    import chemfuse

    entity = entity.strip()
    if not entity:
        raise ValueError(f"No identifier found in query: {original_query!r}")

    logger.debug("nlq.ask: map_identifiers %r", entity)

    # Detect identifier type
    if _is_smiles(entity):
        return chemfuse.map_identifiers(smiles=entity)

    if re.match(r"^CHEMBL\d+$", entity, re.IGNORECASE):
        return chemfuse.map_identifiers(chembl_id=entity.upper())

    if re.match(r"^\d+$", entity):
        return chemfuse.map_identifiers(cid=int(entity))

    # InChIKey: 27-char uppercase with hyphens
    if re.match(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", entity):
        return chemfuse.map_identifiers(inchikey=entity)

    # Try name search and extract InChIKey for mapping
    coll = chemfuse.search(entity, query_type="name")
    if coll and coll[0].inchikey:
        return chemfuse.map_identifiers(inchikey=coll[0].inchikey)

    raise ValueError(f"Cannot determine identifier type for {entity!r}.")


def _handle_compare(entity_a: str, entity_b: str, original_query: str) -> dict[str, Any]:
    """Compare two compounds and return a side-by-side property dict."""
    from chemfuse.compute.druglikeness import check_drug_likeness

    def _resolve_smiles(entity: str) -> str:
        """Return a SMILES string for *entity* (already-SMILES or compound name)."""
        smiles_tokens = _extract_smiles_tokens(entity)
        if smiles_tokens:
            return smiles_tokens[0]
        if _is_smiles(entity):
            return entity
        # Compound name
        import chemfuse
        coll = chemfuse.search(entity, query_type="name")
        if not coll or not coll[0].smiles:
            raise ValueError(f"Cannot resolve SMILES for {entity!r}.")
        return coll[0].smiles

    smiles_a = _resolve_smiles(entity_a)
    smiles_b = _resolve_smiles(entity_b)

    logger.debug("nlq.ask: compare %r and %r", smiles_a, smiles_b)

    result: dict[str, Any] = {
        "compound_a": smiles_a,
        "compound_b": smiles_b,
    }

    # Compute properties for both
    for key, smiles in [("properties_a", smiles_a), ("properties_b", smiles_b)]:
        try:
            from chemfuse.compute.descriptors import compute_descriptors
            result[key] = compute_descriptors(smiles)
        except Exception as exc:  # noqa: BLE001
            logger.debug("nlq compare: descriptor computation failed for %r: %s", smiles, exc)
            result[key] = {}

    # Drug-likeness for both
    for key, smiles, prop_key in [
        ("druglikeness_a", smiles_a, "properties_a"),
        ("druglikeness_b", smiles_b, "properties_b"),
    ]:
        try:
            dl = check_drug_likeness(result.get(prop_key, {}), smiles=smiles)
            result[key] = dl
        except Exception as exc:  # noqa: BLE001
            logger.debug("nlq compare: drug-likeness failed for %r: %s", smiles, exc)
            result[key] = None

    return result
