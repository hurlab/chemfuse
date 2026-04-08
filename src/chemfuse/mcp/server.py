"""ChemFuse MCP server implementation.

Exposes ChemFuse cheminformatics functionality as MCP tools for AI agents.
Supports compound search, ADMET prediction, drug-likeness filters, descriptor
computation, scaffold analysis, and compound comparison.

Guards all mcp imports so the module is importable without the mcp package:
    try:
        from chemfuse.mcp.server import main
    except ImportError:
        ...
"""

from __future__ import annotations

import json
import logging
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional mcp import guard
# ---------------------------------------------------------------------------

try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import TextContent, Tool

    _MCP_AVAILABLE = True
except ImportError:  # pragma: no cover
    _MCP_AVAILABLE = False
    Server = None  # type: ignore[assignment,misc]
    stdio_server = None  # type: ignore[assignment]
    TextContent = None  # type: ignore[assignment]
    Tool = None  # type: ignore[assignment]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _compound_to_dict(compound: Any) -> dict[str, Any]:
    """Convert a Compound model to a JSON-serializable dict."""
    props = compound.properties
    result: dict[str, Any] = {
        "name": compound.name,
        "smiles": compound.smiles,
        "inchikey": compound.inchikey,
        "inchi": compound.inchi,
        "formula": compound.formula,
        "sources": compound.sources,
    }
    if compound.cid is not None:
        result["cid"] = compound.cid
    if compound.chembl_id is not None:
        result["chembl_id"] = compound.chembl_id
    if compound.drugbank_id is not None:
        result["drugbank_id"] = compound.drugbank_id
    if compound.synonyms:
        result["synonyms"] = compound.synonyms[:5]
    if props:
        prop_dict: dict[str, Any] = {}
        if props.molecular_weight is not None:
            prop_dict["molecular_weight"] = props.molecular_weight
        if props.xlogp is not None:
            prop_dict["logp"] = props.xlogp
        if props.tpsa is not None:
            prop_dict["tpsa"] = props.tpsa
        if props.hbd_count is not None:
            prop_dict["hbd"] = props.hbd_count
        if props.hba_count is not None:
            prop_dict["hba"] = props.hba_count
        if props.rotatable_bonds is not None:
            prop_dict["rotatable_bonds"] = props.rotatable_bonds
        if props.heavy_atom_count is not None:
            prop_dict["heavy_atom_count"] = props.heavy_atom_count
        if prop_dict:
            result["properties"] = prop_dict
    return result


def _ok(data: Any) -> list[Any]:
    """Wrap data in a TextContent list."""
    return [TextContent(type="text", text=json.dumps(data, indent=2))]


def _err(message: str) -> list[Any]:
    """Wrap an error message in a TextContent list."""
    return [TextContent(type="text", text=json.dumps({"error": message}))]


# ---------------------------------------------------------------------------
# Tool definitions
# ---------------------------------------------------------------------------

_TOOLS: list[Any] = []

if _MCP_AVAILABLE:
    _TOOLS = [
        Tool(
            name="search_compounds",
            description=(
                "Search for chemical compounds across databases (PubChem, ChEMBL, etc.). "
                "Returns a list of matching compounds with identifiers and properties."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query: compound name, SMILES, CID, formula, or InChI",
                    },
                    "sources": {
                        "type": "array",
                        "items": {"type": "string"},
                        "default": ["pubchem"],
                        "description": "Database sources to search (pubchem, chembl)",
                    },
                    "query_type": {
                        "type": "string",
                        "enum": ["name", "smiles", "cid", "formula", "inchi"],
                        "default": "name",
                        "description": "Type of the query string",
                    },
                    "limit": {
                        "type": "integer",
                        "default": 10,
                        "maximum": 100,
                        "description": "Maximum number of results",
                    },
                },
                "required": ["query"],
            },
        ),
        Tool(
            name="get_compound",
            description=(
                "Retrieve a specific compound by its database identifier. "
                "Use PubChem CIDs (e.g. '2244') for pubchem source."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "identifier": {
                        "type": "string",
                        "description": "Database-specific compound identifier (e.g. PubChem CID)",
                    },
                    "source": {
                        "type": "string",
                        "default": "pubchem",
                        "description": "Source database (pubchem, chembl)",
                    },
                },
                "required": ["identifier"],
            },
        ),
        Tool(
            name="find_similar",
            description=(
                "Find structurally similar compounds using Tanimoto fingerprint similarity. "
                "Returns compounds with similarity scores above the given threshold."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": "SMILES string of the query compound",
                    },
                    "threshold": {
                        "type": "integer",
                        "default": 90,
                        "minimum": 1,
                        "maximum": 100,
                        "description": "Tanimoto similarity threshold (0-100)",
                    },
                    "max_results": {
                        "type": "integer",
                        "default": 10,
                        "maximum": 100,
                        "description": "Maximum number of similar compounds to return",
                    },
                },
                "required": ["smiles"],
            },
        ),
        Tool(
            name="cross_reference",
            description=(
                "Map a compound identifier across multiple databases via UniChem. "
                "Returns a dict of database names to corresponding identifiers."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "identifier": {
                        "type": "string",
                        "description": "The compound identifier value",
                    },
                    "identifier_type": {
                        "type": "string",
                        "enum": ["cid", "chembl_id", "inchikey", "smiles"],
                        "description": "Type of the identifier",
                    },
                },
                "required": ["identifier", "identifier_type"],
            },
        ),
        Tool(
            name="predict_admet",
            description=(
                "Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) "
                "properties for a molecule given its SMILES. "
                "Returns predictions for solubility, BBB permeability, CYP inhibition, hERG "
                "liability, AMES mutagenicity, and an overall drug-likeness score."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": "SMILES string of the molecule",
                    },
                },
                "required": ["smiles"],
            },
        ),
        Tool(
            name="check_drug_likeness",
            description=(
                "Run all standard drug-likeness filters on a molecule: "
                "Lipinski Ro5, Veber, Ghose, Egan, Muegge, PAINS, and QED score."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": "SMILES string of the molecule",
                    },
                },
                "required": ["smiles"],
            },
        ),
        Tool(
            name="compute_descriptors",
            description=(
                "Compute physicochemical molecular descriptors for a SMILES string. "
                "Returns key descriptors: MW, LogP, TPSA, HBD, HBA, rotatable bonds, "
                "ring count, aromatic ring count, QED, and more."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": "SMILES string of the molecule",
                    },
                },
                "required": ["smiles"],
            },
        ),
        Tool(
            name="analyze_scaffolds",
            description=(
                "Analyze Bemis-Murcko scaffolds in a set of compounds. "
                "Returns scaffold SMILES mapped to their frequency count, "
                "sorted by most common scaffold first."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles_list": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of SMILES strings to analyze",
                        "minItems": 1,
                    },
                },
                "required": ["smiles_list"],
            },
        ),
        Tool(
            name="compare_compounds",
            description=(
                "Compare two compounds by their SMILES strings. "
                "Returns Tanimoto similarity, property comparison, and whether "
                "they share the same Murcko scaffold."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "smiles_a": {
                        "type": "string",
                        "description": "SMILES string of the first compound",
                    },
                    "smiles_b": {
                        "type": "string",
                        "description": "SMILES string of the second compound",
                    },
                },
                "required": ["smiles_a", "smiles_b"],
            },
        ),
    ]


# ---------------------------------------------------------------------------
# Handler implementations (pure async functions, testable without Server)
# ---------------------------------------------------------------------------


async def handle_search_compounds(arguments: dict[str, Any]) -> list[Any]:
    """Handle search_compounds tool call."""
    from chemfuse import search_async

    query = arguments["query"]
    sources = arguments.get("sources", ["pubchem"])
    query_type = arguments.get("query_type", "name")
    limit = int(arguments.get("limit", 10))

    try:
        collection = await search_async(query, sources=sources, query_type=query_type, limit=limit)
        results = [_compound_to_dict(c) for c in collection.compounds[:limit]]
        return _ok({"count": len(results), "compounds": results})
    except Exception as exc:  # noqa: BLE001
        logger.exception("search_compounds failed for query=%r", query)
        return _err(f"Search failed: {exc}")


async def handle_get_compound(arguments: dict[str, Any]) -> list[Any]:
    """Handle get_compound tool call."""
    from chemfuse import get_async

    identifier = arguments["identifier"]
    source = arguments.get("source", "pubchem")

    try:
        compound = await get_async(identifier, source=source)
        if compound is None:
            return _err(f"Compound '{identifier}' not found in {source}")
        return _ok(_compound_to_dict(compound))
    except Exception as exc:  # noqa: BLE001
        logger.exception("get_compound failed for identifier=%r", identifier)
        return _err(f"Retrieval failed: {exc}")


async def handle_find_similar(arguments: dict[str, Any]) -> list[Any]:
    """Handle find_similar tool call."""
    from chemfuse import find_similar_async

    smiles = arguments["smiles"]
    threshold = int(arguments.get("threshold", 90))
    max_results = int(arguments.get("max_results", 10))

    try:
        collection = await find_similar_async(smiles, threshold=threshold, max_results=max_results)
        results = [_compound_to_dict(c) for c in collection.compounds[:max_results]]
        return _ok({"count": len(results), "query_smiles": smiles, "compounds": results})
    except Exception as exc:  # noqa: BLE001
        logger.exception("find_similar failed for smiles=%r", smiles)
        return _err(f"Similarity search failed: {exc}")


async def handle_cross_reference(arguments: dict[str, Any]) -> list[Any]:
    """Handle cross_reference tool call."""
    from chemfuse import map_identifiers_async

    identifier = arguments["identifier"]
    identifier_type = arguments["identifier_type"]

    kwargs: dict[str, Any] = {}
    if identifier_type == "cid":
        try:
            kwargs["cid"] = int(identifier)
        except ValueError:
            return _err(f"CID must be an integer, got: {identifier!r}")
    elif identifier_type == "chembl_id":
        kwargs["chembl_id"] = identifier
    elif identifier_type == "inchikey":
        kwargs["inchikey"] = identifier
    elif identifier_type == "smiles":
        kwargs["smiles"] = identifier
    else:
        return _err(f"Unknown identifier_type: {identifier_type!r}")

    try:
        mappings = await map_identifiers_async(**kwargs)
        return _ok({"identifier": identifier, "identifier_type": identifier_type, "mappings": mappings})
    except Exception as exc:  # noqa: BLE001
        logger.exception("cross_reference failed for %s=%r", identifier_type, identifier)
        return _err(f"Cross-reference failed: {exc}")


def handle_predict_admet(arguments: dict[str, Any]) -> list[Any]:
    """Handle predict_admet tool call (synchronous — uses RDKit locally)."""
    from chemfuse.compute.admet import predict_admet
    from chemfuse.core.exceptions import ChemFuseError

    smiles = arguments["smiles"]

    try:
        profile = predict_admet(smiles)
        preds_dict: dict[str, Any] = {}
        for name, pred in profile.predictions.items():
            preds_dict[name] = {
                "value": pred.value,
                "category": pred.category,
                "method": pred.method,
            }
            if pred.unit:
                preds_dict[name]["unit"] = pred.unit
            if pred.confidence is not None:
                preds_dict[name]["confidence"] = pred.confidence
        return _ok({
            "smiles": profile.smiles,
            "overall_score": profile.overall_score,
            "predictions": preds_dict,
        })
    except ChemFuseError as exc:
        return _err(f"Invalid SMILES or prediction error: {exc}")
    except Exception as exc:  # noqa: BLE001
        logger.exception("predict_admet failed for smiles=%r", smiles)
        return _err(f"ADMET prediction failed: {exc}")


def handle_check_drug_likeness(arguments: dict[str, Any]) -> list[Any]:
    """Handle check_drug_likeness tool call (synchronous — uses RDKit locally)."""
    from chemfuse.compute.druglikeness import check_drug_likeness

    smiles = arguments["smiles"]

    try:
        # Pass empty properties dict; check_drug_likeness fills them from SMILES via RDKit
        dl = check_drug_likeness(properties={}, smiles=smiles)
        result: dict[str, Any] = {
            "smiles": smiles,
            "overall_pass": dl.overall_pass,
            "filters": {
                "lipinski": {
                    "pass": dl.lipinski.pass_filter,
                    "violations": dl.lipinski.violations,
                    "details": dl.lipinski.details,
                },
                "veber": {
                    "pass": dl.veber.pass_filter,
                    "violations": dl.veber.violations,
                    "details": dl.veber.details,
                },
                "ghose": {
                    "pass": dl.ghose.pass_filter,
                    "violations": dl.ghose.violations,
                    "details": dl.ghose.details,
                },
                "egan": {
                    "pass": dl.egan.pass_filter,
                    "violations": dl.egan.violations,
                    "details": dl.egan.details,
                },
                "muegge": {
                    "pass": dl.muegge.pass_filter,
                    "violations": dl.muegge.violations,
                    "details": dl.muegge.details,
                },
            },
        }
        if dl.pains is not None:
            result["filters"]["pains"] = {
                "pass": dl.pains.pass_filter,
                "violations": dl.pains.violations,
            }
        if dl.qed is not None:
            result["qed"] = dl.qed
        return _ok(result)
    except Exception as exc:  # noqa: BLE001
        logger.exception("check_drug_likeness failed for smiles=%r", smiles)
        return _err(f"Drug-likeness check failed: {exc}")


def handle_compute_descriptors(arguments: dict[str, Any]) -> list[Any]:
    """Handle compute_descriptors tool call (synchronous — uses RDKit locally)."""
    from chemfuse.compute.descriptors import compute_descriptors

    smiles = arguments["smiles"]

    # Keys extracted for the summary response sent to AI agents.
    # The full RDKit descriptor set (200+) is available but too verbose for
    # agent consumption; select the most clinically relevant ones.
    _KEY_DESCRIPTORS = {
        "MolWt", "ExactMolWt", "MolLogP", "TPSA",
        "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
        "NumAromaticRings", "NumRings", "HeavyAtomCount",
        "FractionCSP3", "qed",
        "BertzCT", "BalabanJ",
        "NumSaturatedRings", "NumAliphaticRings",
        "MolMR", "LabuteASA",
        "MaxAbsPartialCharge", "MinAbsPartialCharge",
    }

    try:
        all_descriptors = compute_descriptors(smiles)
        if not all_descriptors:
            return _err(f"Could not compute descriptors for SMILES: {smiles!r}")

        # Return selected key descriptors plus total count
        summary = {k: v for k, v in all_descriptors.items() if k in _KEY_DESCRIPTORS}
        return _ok({
            "smiles": smiles,
            "total_descriptors": len(all_descriptors),
            "descriptors": summary,
        })
    except ImportError:
        return _err("RDKit is required. Install with: pip install chemfuse[rdkit]")
    except Exception as exc:  # noqa: BLE001
        logger.exception("compute_descriptors failed for smiles=%r", smiles)
        return _err(f"Descriptor computation failed: {exc}")


def handle_analyze_scaffolds(arguments: dict[str, Any]) -> list[Any]:
    """Handle analyze_scaffolds tool call (synchronous — uses RDKit locally)."""
    from chemfuse.analyze.scaffolds import scaffold_frequency
    from chemfuse.core.exceptions import OptionalDependencyError

    smiles_list: list[str] = arguments["smiles_list"]
    if not smiles_list:
        return _err("smiles_list must not be empty")

    try:
        freq = scaffold_frequency(smiles_list)
        scaffolds = [
            {"scaffold_smiles": smi, "count": count, "frequency": round(count / len(smiles_list), 3)}
            for smi, count in freq.items()
        ]
        no_scaffold = len(smiles_list) - sum(s["count"] for s in scaffolds)
        return _ok({
            "total_compounds": len(smiles_list),
            "unique_scaffolds": len(scaffolds),
            "no_scaffold": no_scaffold,
            "scaffolds": scaffolds,
        })
    except OptionalDependencyError:
        return _err("RDKit is required. Install with: pip install chemfuse[rdkit]")
    except Exception as exc:  # noqa: BLE001
        logger.exception("analyze_scaffolds failed")
        return _err(f"Scaffold analysis failed: {exc}")


def handle_compare_compounds(arguments: dict[str, Any]) -> list[Any]:
    """Handle compare_compounds tool call (synchronous — uses RDKit locally)."""
    smiles_a = arguments["smiles_a"]
    smiles_b = arguments["smiles_b"]

    try:
        from chemfuse.compute.fingerprints import tanimoto_similarity
        tanimoto = tanimoto_similarity(smiles_a, smiles_b)
    except ImportError:
        return _err("RDKit is required. Install with: pip install chemfuse[rdkit]")
    except ValueError as exc:
        return _err(f"Invalid SMILES: {exc}")
    except Exception as exc:  # noqa: BLE001
        return _err(f"Similarity computation failed: {exc}")

    # Property comparison via descriptors
    property_comparison: dict[str, Any] = {}
    try:
        from chemfuse.compute.descriptors import compute_descriptors

        _COMPARE_KEYS = ["MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds"]
        desc_a = compute_descriptors(smiles_a)
        desc_b = compute_descriptors(smiles_b)
        for key in _COMPARE_KEYS:
            val_a = desc_a.get(key)
            val_b = desc_b.get(key)
            if val_a is not None and val_b is not None:
                property_comparison[key] = {
                    "compound_a": round(val_a, 3),
                    "compound_b": round(val_b, 3),
                    "delta": round(val_b - val_a, 3),
                }
    except Exception:  # noqa: BLE001
        pass  # property comparison is best-effort

    # Scaffold comparison
    shared_scaffold: str | None = None
    try:
        from chemfuse.analyze.scaffolds import murcko_scaffold

        scaffold_a = murcko_scaffold(smiles_a)
        scaffold_b = murcko_scaffold(smiles_b)
        if scaffold_a and scaffold_b and scaffold_a == scaffold_b:
            shared_scaffold = scaffold_a
    except Exception:  # noqa: BLE001
        pass  # scaffold comparison is best-effort

    return _ok({
        "smiles_a": smiles_a,
        "smiles_b": smiles_b,
        "tanimoto_similarity": round(tanimoto, 4),
        "similarity_label": (
            "high" if tanimoto >= 0.7 else "medium" if tanimoto >= 0.4 else "low"
        ),
        "shared_scaffold": shared_scaffold,
        "property_comparison": property_comparison,
    })


# ---------------------------------------------------------------------------
# Server setup
# ---------------------------------------------------------------------------

if _MCP_AVAILABLE:
    # @MX:ANCHOR: [AUTO] Public MCP server instance — entry point for all tool calls
    # @MX:REASON: Central server object referenced by __main__ and list_tools/call_tool handlers
    app = Server("chemfuse")

    @app.list_tools()
    async def list_tools() -> list[Any]:
        """Return all available ChemFuse tools."""
        return _TOOLS

    @app.call_tool()
    async def call_tool(name: str, arguments: dict[str, Any]) -> list[Any]:
        """Dispatch tool calls to the appropriate handler."""
        # Async handlers (network I/O via ChemFuse API clients)
        if name == "search_compounds":
            return await handle_search_compounds(arguments)
        if name == "get_compound":
            return await handle_get_compound(arguments)
        if name == "find_similar":
            return await handle_find_similar(arguments)
        if name == "cross_reference":
            return await handle_cross_reference(arguments)

        # Synchronous handlers (local RDKit computation)
        if name == "predict_admet":
            return handle_predict_admet(arguments)
        if name == "check_drug_likeness":
            return handle_check_drug_likeness(arguments)
        if name == "compute_descriptors":
            return handle_compute_descriptors(arguments)
        if name == "analyze_scaffolds":
            return handle_analyze_scaffolds(arguments)
        if name == "compare_compounds":
            return handle_compare_compounds(arguments)

        return _err(f"Unknown tool: {name!r}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


async def main() -> None:
    """Run the ChemFuse MCP server over stdio."""
    if not _MCP_AVAILABLE:
        raise ImportError(
            "The 'mcp' package is required to run the MCP server. "
            "Install it with: pip install chemfuse[mcp]"
        )
    async with stdio_server() as (read_stream, write_stream):
        await app.run(read_stream, write_stream, app.create_initialization_options())
