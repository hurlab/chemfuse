"""Search page - multi-source compound search with molecule grid."""

from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st

from chemfuse.web._utils import _run_async

# Maximum number of compounds to store in session to prevent unbounded growth
MAX_SESSION_COMPOUNDS = 100

# Search type labels map to actual query_type values supported by adapters
_SEARCH_TYPE_OPTIONS = ["name", "smiles", "cid", "formula", "inchi"]
_SEARCH_TYPE_LABELS = {
    "name": "Name",
    "smiles": "SMILES",
    "cid": "CID",
    "formula": "Formula",
    "inchi": "InChI",
}


def render() -> None:
    """Render the Search page."""
    st.header("Compound Search")
    st.caption("Search across PubChem, ChEMBL, and more.")

    # -- Search controls --
    col_query, col_type = st.columns([3, 1])

    with col_query:
        query = st.text_input(
            "Query",
            placeholder="e.g., aspirin  |  CC(=O)Oc1ccccc1C(=O)O  |  2244",
            key="search_query_input",
        )

    with col_type:
        search_type = st.selectbox(
            "Search type",
            _SEARCH_TYPE_OPTIONS,
            format_func=lambda k: _SEARCH_TYPE_LABELS.get(k, k),
            key="search_type_select",
        )

    # Source selection
    st.markdown("**Sources**")
    col_pc, col_chembl = st.columns(2)
    with col_pc:
        use_pubchem = st.checkbox("PubChem", value=True, key="src_pubchem")
    with col_chembl:
        use_chembl = st.checkbox("ChEMBL", value=False, key="src_chembl")

    search_clicked = st.button("Search", type="primary", use_container_width=True, key="search_btn")

    # Execute search
    if search_clicked and query.strip():
        sources = []
        if use_pubchem:
            sources.append("pubchem")
        if use_chembl:
            sources.append("chembl")
        if not sources:
            st.warning("Select at least one source.")
            return
        _execute_search(query.strip(), search_type, sources)

    # Display persisted results
    results = st.session_state.get("search_results", [])
    if results:
        _display_results(results)


def _execute_search(query: str, search_type: str, sources: list[str]) -> None:
    """Execute search and store results in session_state."""
    with st.spinner(f"Searching for '{query}'..."):
        # Convert list to tuple for @st.cache_data hashability (C7)
        results, warnings = _fetch_results(query, search_type, tuple(sources))

    # Display any warnings/errors from the cached function (C6)
    for w in warnings:
        st.warning(w)

    if not results:
        st.warning("No results found.")
        st.session_state["search_results"] = []
        return

    st.session_state["search_results"] = results
    # Offer to add results to session compounds for comparison
    compounds = st.session_state.get("session_compounds", [])
    existing_smiles = {c.get("smiles", "") for c in compounds}
    new_compounds = [r for r in results if r.get("smiles", "") not in existing_smiles]
    if new_compounds:
        compounds.extend(new_compounds[:20])
        # Cap session_compounds to prevent unbounded growth (C12)
        if len(compounds) > MAX_SESSION_COMPOUNDS:
            compounds = compounds[-MAX_SESSION_COMPOUNDS:]
        st.session_state["session_compounds"] = compounds


@st.cache_data(ttl=300, show_spinner=False)
def _fetch_results(
    query: str,
    search_type: str,
    sources: tuple[str, ...],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Fetch compound search results (cached, TTL=5 min).

    No st.* calls are made inside this function to comply with Streamlit's
    caching requirements (C6).

    Args:
        query: The search query string.
        search_type: One of the supported query types (name, smiles, cid, ...).
        sources: Tuple of source names (e.g. ('pubchem', 'chembl')).

    Returns:
        A tuple of (results_list, warnings_list).
    """
    warnings: list[str] = []

    async def _run() -> list[dict[str, Any]]:
        collected: list[dict[str, Any]] = []
        for source_name in sources:
            try:
                if source_name == "pubchem":
                    from chemfuse.sources.pubchem import PubChemSource
                    src = PubChemSource()
                    async with src:
                        compounds = await src.search(query, query_type=search_type)
                    collected.extend(_compound_to_dict(c) for c in compounds)
                elif source_name == "chembl":
                    from chemfuse.sources.chembl import ChEMBLSource
                    src = ChEMBLSource()
                    async with src:
                        compounds = await src.search(query, query_type=search_type)
                    collected.extend(_compound_to_dict(c) for c in compounds)
            except Exception as exc:
                # Collect warnings instead of calling st.warning (C6)
                warnings.append(f"Source '{source_name}' error: {exc}")
        return collected

    try:
        results = _run_async(_run())
    except Exception as exc:
        warnings.append(f"Search failed: {exc}")
        results = []

    return results, warnings


def _compound_to_dict(compound: Any) -> dict[str, Any]:
    """Convert a Compound model to a plain dictionary."""
    if hasattr(compound, "to_dict"):
        return compound.to_dict()
    if hasattr(compound, "model_dump"):
        return compound.model_dump()
    return dict(compound) if compound else {}


def _display_results(results: list[dict[str, Any]]) -> None:
    """Display search results as a molecule grid with export."""
    st.subheader(f"Results ({len(results)} compounds)")

    col_view, col_export = st.columns([2, 1])
    with col_view:
        view_mode = st.radio("View", ["Grid", "Table"], horizontal=True, key="search_view_mode")
    with col_export:
        if st.button("Export CSV", key="export_csv_btn"):
            _export_csv(results)

    if view_mode == "Grid":
        from chemfuse.web.components.mol_grid import render_mol_grid
        render_mol_grid(results, on_click_session_key="selected_compound")
    else:
        _display_table(results)


def _display_table(results: list[dict[str, Any]]) -> None:
    """Display results as a sortable dataframe."""
    df = pd.DataFrame(results)
    preferred_cols = ["name", "cid", "chembl_id", "formula", "smiles",
                      "molecular_weight", "xlogp", "tpsa"]
    display_cols = [c for c in preferred_cols if c in df.columns]
    if not display_cols:
        display_cols = list(df.columns)
    st.dataframe(df[display_cols], use_container_width=True, hide_index=True)


def _export_csv(results: list[dict[str, Any]]) -> None:
    """Trigger a CSV download of search results."""
    df = pd.DataFrame(results)
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)
    st.download_button(
        label="Download CSV",
        data=buf,
        file_name="chemfuse_search_results.csv",
        mime="text/csv",
        key="download_csv_btn",
    )
