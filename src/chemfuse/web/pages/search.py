"""Search page - multi-source compound search with molecule grid."""

from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st


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
            ["name", "smiles", "similarity", "substructure"],
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
        results = _fetch_results(query, search_type, sources)

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
        st.session_state["session_compounds"] = compounds


@st.cache_data(ttl=300, show_spinner=False)
def _fetch_results(query: str, search_type: str, sources: tuple[str, ...]) -> list[dict[str, Any]]:
    """Fetch compound search results (cached, TTL=5 min).

    Uses synchronous HTTP via httpx to avoid Streamlit's async limitations.
    Falls back gracefully if sources are unavailable.
    """
    import asyncio

    results: list[dict[str, Any]] = []

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
                st.warning(f"Source '{source_name}' error: {exc}")
        return collected

    try:
        loop = asyncio.new_event_loop()
        results = loop.run_until_complete(_run())
        loop.close()
    except Exception as exc:
        st.error(f"Search failed: {exc}")

    return results


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
