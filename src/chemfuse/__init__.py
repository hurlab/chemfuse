"""ChemFuse: Multi-database cheminformatics suite.

Public API for searching chemical databases, retrieving compound data,
and performing cross-database analysis.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any

from chemfuse._version import __version__
from chemfuse.core._async import run_async as _run_sync
from chemfuse.core.exceptions import NotFoundError, SourceError
from chemfuse.models.bioactivity import BindingMeasurement, Bioactivity
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties
from chemfuse.models.patent import Patent
from chemfuse.models.prediction import ADMETPrediction, ADMETProfile, DrugLikeness, FilterResult
from chemfuse.models.target import TargetAssociation
from chemfuse.nlq import ask
from chemfuse.sources import registry

logger = logging.getLogger(__name__)


def _try_standardize_for_dedup(compounds: list[Compound]) -> None:
    """Attempt to standardize SMILES and fill missing InChIKeys for better dedup.

    Operates in-place. Failures are silently skipped (best-effort).
    """
    try:
        from rdkit import Chem  # noqa: I001
        from chemfuse.compute.standardization import standardize_mol
    except ImportError:
        return  # RDKit not available, skip standardization

    for compound in compounds:
        if compound.smiles and not compound.inchikey:
            try:
                std_smiles = standardize_mol(compound.smiles)
                if std_smiles:
                    mol = Chem.MolFromSmiles(std_smiles)
                    if mol:
                        from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi  # noqa: I001
                        inchi = MolToInchi(mol)
                        if inchi:
                            compound.inchikey = InchiToInchiKey(inchi)
            except Exception:  # noqa: BLE001
                pass  # Best-effort: skip this compound


def _merge_by_inchikey(compounds: list[Compound]) -> list[Compound]:
    """Deduplicate compounds by merging those with the same InChIKey.

    Applies tautomer canonicalization before InChIKey comparison to prevent
    duplicate entries from different tautomeric forms of the same compound.

    Args:
        compounds: Flat list of compounds from all sources.

    Returns:
        Deduplicated list with merged Compound objects.
    """
    # Try to standardize SMILES and recompute InChIKey for better matching
    _try_standardize_for_dedup(compounds)

    # Group by InChIKey; compounds without InChIKey are kept as-is
    keyed: dict[str, Compound] = {}
    no_key: list[Compound] = []

    for compound in compounds:
        if compound.inchikey:
            key = compound.inchikey
            if key in keyed:
                try:
                    keyed[key] = keyed[key].merge(compound)
                except ValueError:
                    # Different molecule; keep both (should not happen for same InChIKey)
                    no_key.append(compound)
            else:
                keyed[key] = compound
        else:
            no_key.append(compound)

    return list(keyed.values()) + no_key


async def search_async(
    query: str,
    sources: list[str] | None = None,
    query_type: str = "name",
    limit: int = 100,
    progress_callback: Any = None,
) -> CompoundCollection:
    """Search for compounds across one or more chemical databases (async).

    Searches all specified sources in parallel using asyncio.gather().
    Results from multiple sources are merged by InChIKey.

    Args:
        query: Search query (name, SMILES, CID, formula, InChI).
        sources: List of source names to query. Defaults to ['pubchem'].
        query_type: Query type: 'name', 'smiles', 'cid', 'formula', 'inchi'.
        limit: Maximum results per source.
        progress_callback: Optional callable(source_name) called before each source.

    Returns:
        CompoundCollection with merged, deduplicated results from all sources.
    """
    if sources is None:
        sources = ["pubchem"]

    warnings: list[str] = []

    async def _search_one(source_name: str) -> list[Compound]:
        adapter = registry.get(source_name)
        if progress_callback is not None:
            progress_callback(source_name)
        return await adapter.search(query, query_type=query_type)

    # Run all sources in parallel
    results = await asyncio.gather(
        *[_search_one(s) for s in sources],
        return_exceptions=True,
    )

    all_compounds: list[Compound] = []
    successful_sources: list[str] = []

    for source_name, result in zip(sources, results, strict=True):
        if isinstance(result, Exception):
            warnings.append(f"{source_name}: {result}")
            logger.warning("search: source %s failed: %s", source_name, result)
        else:
            all_compounds.extend(result[:limit])
            successful_sources.append(source_name)

    # Merge compounds by InChIKey
    merged = _merge_by_inchikey(all_compounds)

    collection = CompoundCollection(
        compounds=merged,
        query=query,
        sources=successful_sources,
    )
    # Attach warnings to collection for inspection by caller
    if warnings:
        collection.warnings = warnings

    return collection


def search(
    query: str,
    sources: list[str] | None = None,
    query_type: str = "name",
    limit: int = 100,
) -> CompoundCollection:
    """Search for compounds across one or more chemical databases (sync).

    Args:
        query: Search query (name, SMILES, CID, formula, InChI).
        sources: List of source names to query. Defaults to ['pubchem'].
        query_type: Query type: 'name', 'smiles', 'cid', 'formula', 'inchi'.
        limit: Maximum results per source.

    Returns:
        CompoundCollection with merged results from all sources.
    """
    return _run_sync(search_async(query, sources=sources, query_type=query_type, limit=limit))


async def get_async(
    identifier: str,
    source: str = "pubchem",
) -> Compound | None:
    """Retrieve a specific compound by identifier (async).

    Args:
        identifier: Database-specific compound identifier (e.g., PubChem CID).
        source: Source database name.

    Returns:
        Compound object, or None if not found.
    """
    adapter = registry.get(source)
    return await adapter.get_by_id(identifier)


def get(
    identifier: str,
    source: str = "pubchem",
) -> Compound | None:
    """Retrieve a specific compound by identifier (sync).

    Args:
        identifier: Database-specific compound identifier (e.g., PubChem CID).
        source: Source database name.

    Returns:
        Compound object, or None if not found.
    """
    return _run_sync(get_async(identifier, source=source))


async def find_similar_async(
    smiles: str,
    threshold: int = 90,
    max_results: int = 100,
    sources: list[str] | None = None,
) -> CompoundCollection:
    """Find compounds similar to a given SMILES string (async).

    Args:
        smiles: SMILES string for the query compound.
        threshold: Tanimoto similarity threshold (0-100).
        max_results: Maximum number of results per source.
        sources: List of source names. Defaults to ['pubchem'].

    Returns:
        CompoundCollection with similar compounds.
    """
    if sources is None:
        sources = ["pubchem"]

    # Filter to only sources that support get_similarity
    capable_sources = []
    for source_name in sources:
        adapter = registry.get(source_name)
        if getattr(adapter, "get_similarity", None) is not None:
            capable_sources.append(source_name)

    async def _search_similar_one(source_name: str) -> list[Compound]:
        adapter = registry.get(source_name)
        get_sim = adapter.get_similarity
        return await get_sim(smiles, threshold=threshold, max_results=max_results)

    results = await asyncio.gather(
        *[_search_similar_one(s) for s in capable_sources],
        return_exceptions=True,
    )

    all_compounds: list[Compound] = []
    source_names: list[str] = []

    for source_name, result in zip(capable_sources, results, strict=True):
        if isinstance(result, (SourceError, NotFoundError, Exception)):
            logger.debug("find_similar_async: source %s failed: %s", source_name, result)
        else:
            all_compounds.extend(result)
            if source_name not in source_names:
                source_names.append(source_name)

    return CompoundCollection(
        compounds=all_compounds,
        query=smiles,
        sources=source_names,
    )


def find_similar(
    smiles: str,
    threshold: int = 90,
    max_results: int = 100,
    sources: list[str] | None = None,
) -> CompoundCollection:
    """Find compounds similar to a given SMILES string (sync).

    Args:
        smiles: SMILES string for the query compound.
        threshold: Tanimoto similarity threshold (0-100).
        max_results: Maximum number of results per source.
        sources: List of source names. Defaults to ['pubchem'].

    Returns:
        CompoundCollection with similar compounds.
    """
    return _run_sync(
        find_similar_async(
            smiles,
            threshold=threshold,
            max_results=max_results,
            sources=sources,
        )
    )


async def cross_reference_async(
    compound: Compound,
    target_sources: list[str] | None = None,
) -> Compound:
    """Look up a compound across multiple databases and merge results (async).

    Uses InChIKey or SMILES as the cross-reference key.

    Args:
        compound: Source compound with at least a SMILES or InChIKey.
        target_sources: List of source names to cross-reference.

    Returns:
        Merged Compound with data from all sources.
    """
    if target_sources is None:
        target_sources = ["pubchem"]

    merged = compound

    for source_name in target_sources:
        adapter = registry.get(source_name)
        if compound.inchikey:
            query = compound.inchikey
            query_type = "inchikey"
        elif compound.smiles:
            query = compound.smiles
            query_type = "smiles"
        else:
            query = compound.name
            query_type = "name"
        if not query:
            continue
        try:
            results = await adapter.search(query, query_type=query_type)
            if results:
                merged = merged.merge(results[0])
        except Exception as e:  # noqa: BLE001
            logger.debug("cross_reference_async: source %s failed: %s", source_name, e)

    return merged


def cross_reference(
    compound: Compound,
    target_sources: list[str] | None = None,
) -> Compound:
    """Look up a compound across multiple databases and merge results (sync).

    Args:
        compound: Source compound with at least a SMILES or InChIKey.
        target_sources: List of source names to cross-reference.

    Returns:
        Merged Compound with data from all sources.
    """
    return _run_sync(cross_reference_async(compound, target_sources=target_sources))


async def map_identifiers_async(
    cid: int | None = None,
    chembl_id: str | None = None,
    inchikey: str | None = None,
    smiles: str | None = None,
) -> dict[str, str]:
    """Map a compound identifier to other database IDs using UniChem (async).

    At least one of cid, chembl_id, inchikey, or smiles must be provided.

    Args:
        cid: PubChem CID.
        chembl_id: ChEMBL ID (e.g., "CHEMBL25").
        inchikey: Standard InChIKey (27 chars).
        smiles: SMILES string (resolved to InChIKey via PubChem first).

    Returns:
        Dictionary mapping database names to identifiers.
        Returns dict with only the input entry if no mappings found.
    """
    from chemfuse.sources import registry

    unichem = registry.get("unichem")
    result: dict[str, str] = {}

    if inchikey:
        result = await unichem.cross_reference(inchikey)
    elif cid is not None:
        result = await unichem.map_identifiers(str(cid), "pubchem")
        if not result:
            result["pubchem"] = str(cid)
    elif chembl_id:
        result = await unichem.map_identifiers(chembl_id, "chembl")
        if not result:
            result["chembl"] = chembl_id
    elif smiles:
        # Resolve SMILES to InChIKey via PubChem
        try:
            pubchem = registry.get("pubchem")
            compounds = await pubchem.search(smiles, query_type="smiles")
            if compounds and compounds[0].inchikey:
                result = await unichem.cross_reference(compounds[0].inchikey)
        except Exception as exc:
            logger.debug("map_identifiers: SMILES lookup failed: %s", exc)

    if not result:
        # Return at minimum the input
        if cid is not None:
            result = {"pubchem": str(cid)}
        elif chembl_id:
            result = {"chembl": chembl_id}
        elif inchikey:
            result = {"inchikey": inchikey}
        elif smiles:
            result = {"smiles": smiles}

    return result


def map_identifiers(
    cid: int | None = None,
    chembl_id: str | None = None,
    inchikey: str | None = None,
    smiles: str | None = None,
) -> dict[str, str]:
    """Map a compound identifier to other database IDs using UniChem (sync).

    Args:
        cid: PubChem CID.
        chembl_id: ChEMBL ID.
        inchikey: Standard InChIKey.
        smiles: SMILES string.

    Returns:
        Dictionary mapping database names to identifiers.
    """
    return _run_sync(
        map_identifiers_async(cid=cid, chembl_id=chembl_id, inchikey=inchikey, smiles=smiles)
    )


def batch_search(
    input_path: str,
    sources: list[str] | None = None,
    concurrency: int = 5,
    **kwargs: Any,
) -> tuple[CompoundCollection, list[dict[str, Any]]]:
    """Search for compounds from a file in batch (sync).

    Args:
        input_path: Path to CSV or TXT file with queries.
        sources: List of source names. Defaults to ['pubchem'].
        concurrency: Maximum concurrent requests.
        **kwargs: Additional keyword arguments passed to batch_search_async.

    Returns:
        Tuple of (CompoundCollection, list of error dicts).
    """
    from chemfuse.core.batch import batch_search as _batch_search

    return _batch_search(
        input_path,
        sources=sources,
        concurrency=concurrency,
        **kwargs,
    )


# Register pandas accessor when pandas is available
try:
    import chemfuse.pandas_ext  # noqa: F401
except ImportError:
    pass


__all__ = [
    "__version__",
    "ask",
    "search",
    "search_async",
    "get",
    "get_async",
    "find_similar",
    "find_similar_async",
    "cross_reference",
    "cross_reference_async",
    "map_identifiers",
    "map_identifiers_async",
    "batch_search",
    "Compound",
    "CompoundProperties",
    "CompoundCollection",
    "Bioactivity",
    "BindingMeasurement",
    "TargetAssociation",
    "Patent",
    "FilterResult",
    "DrugLikeness",
    "ADMETProfile",
    "ADMETPrediction",
    "NotFoundError",
    "SourceError",
]
