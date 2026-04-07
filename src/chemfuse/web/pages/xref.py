"""Cross-Reference page - identifier mapping across databases."""

from __future__ import annotations

import re
from html import escape
from typing import Any

import pandas as pd
import streamlit as st

from chemfuse.web._utils import _run_async

# Database URL templates for clickable links
_DB_URLS: dict[str, str] = {
    "pubchem": "https://pubchem.ncbi.nlm.nih.gov/compound/{id}",
    "chembl": "https://www.ebi.ac.uk/chembl/compound_report_card/{id}/",
    "bindingdb": "https://www.bindingdb.org/rwd/bind/chemsearch/rwd/AdvanceSearch.jsp?target={id}",
    "drugbank": "https://go.drugbank.com/drugs/{id}",
    "chebi": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{id}",
    "stdinchikey": "https://www.chemspider.com/InChIKey/{id}",
    "hmdb": "https://hmdb.ca/metabolites/{id}",
}

_DB_LABELS: dict[str, str] = {
    "pubchem": "PubChem CID",
    "chembl": "ChEMBL ID",
    "bindingdb": "BindingDB",
    "drugbank": "DrugBank",
    "chebi": "ChEBI",
    "stdinchikey": "InChIKey",
    "hmdb": "HMDB",
    "cas": "CAS Number",
}


def render() -> None:
    """Render the Cross-Reference page."""
    st.header("Cross-Reference")
    st.caption("Map a compound identifier across all supported databases.")

    # Input
    col_id, col_btn = st.columns([3, 1])
    with col_id:
        identifier = st.text_input(
            "Enter identifier",
            key="xref_input",
            placeholder="CID, ChEMBL ID, SMILES, InChIKey, or compound name",
        )
    with col_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        lookup_clicked = st.button("Look up", type="primary", key="xref_lookup_btn")

    if lookup_clicked and identifier.strip():
        with st.spinner("Resolving cross-references..."):
            lookup_return = _lookup_xref(identifier.strip())
        # Handle both tuple return (result, error) and None (e.g. from mocks)
        if isinstance(lookup_return, tuple):
            result, error = lookup_return
        else:
            result, error = lookup_return, None
        if error:
            st.error(error)
        if result:
            st.session_state["xref_result"] = {
                "query": identifier.strip(),
                "mappings": result,
            }
        else:
            st.warning("No cross-references found for this identifier.")
            st.session_state["xref_result"] = None

    # Display result
    xref_data = st.session_state.get("xref_result")
    if xref_data:
        # Support both new format {"query": ..., "mappings": ...}
        # and legacy format {db_key: db_id, ...}
        if "mappings" in xref_data:
            stored_query = xref_data.get("query", identifier.strip())
            _display_xref_table(stored_query, xref_data["mappings"])
        else:
            # Legacy format: the dict itself is the mappings
            _display_xref_table(identifier.strip(), xref_data)


@st.cache_data(ttl=600, show_spinner=False)
def _lookup_xref(identifier: str) -> tuple[dict[str, str] | None, str | None]:
    """Look up cross-references for an identifier via UniChem.

    Returns:
        A tuple of (result_dict_or_None, error_message_or_None).
    """
    async def _run() -> dict[str, str] | None:
        from chemfuse.sources import registry
        src = registry.get("unichem")
        if identifier.isdigit():
            return await src.map_identifiers(identifier, "pubchem")
        elif identifier.upper().startswith("CHEMBL"):
            return await src.map_identifiers(identifier, "chembl")
        elif re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', identifier):
            return await src.cross_reference(identifier)
        else:
            return await _resolve_via_pubchem(identifier, src)

    try:
        result = _run_async(_run())
        return result, None
    except Exception as exc:
        return None, f"Cross-reference lookup failed: {exc}"


async def _resolve_via_pubchem(identifier: str, unichem_src: Any) -> dict[str, str] | None:
    """Resolve identifier via PubChem then get cross-references."""
    try:
        from chemfuse.sources import registry
        pubchem_src = registry.get("pubchem")
        compounds = await pubchem_src.search(identifier, query_type="name")
        if not compounds:
            compounds = await pubchem_src.search(identifier, query_type="smiles")
        if compounds:
            compound = compounds[0]
            cid = compound.cid
            if cid:
                mappings = await unichem_src.map_identifiers(str(cid), "pubchem")
                if mappings:
                    mappings["pubchem"] = str(cid)
                return mappings
    except Exception:
        pass
    return None


def _display_xref_table(query: str, mappings: dict[str, str]) -> None:
    """Display cross-reference table with clickable links."""
    st.subheader(f"Cross-references for: {query}")

    rows: list[dict[str, str]] = []
    for db_key, db_id in mappings.items():
        if not db_id:
            continue
        label = _DB_LABELS.get(db_key, db_key.capitalize())
        url_template = _DB_URLS.get(db_key)
        if url_template:
            # HTML-escape all dynamic values before embedding in HTML
            safe_id = escape(str(db_id))
            link = url_template.format(id=db_id)
            safe_url = escape(link)
            link_html = f'<a href="{safe_url}" target="_blank">{safe_id}</a>'
        else:
            link_html = escape(str(db_id))
        rows.append({
            "Database": label,
            "Identifier": db_id,
            "Link": link_html,
        })

    if not rows:
        st.info("No mappings found.")
        return

    df = pd.DataFrame(rows)

    # Display plain table
    st.dataframe(
        df[["Database", "Identifier"]],
        use_container_width=True,
        hide_index=True,
    )

    # Display clickable links via HTML
    st.markdown("**Database Links**")
    for row in rows:
        if "Link" in row and "<a href" in row["Link"]:
            st.markdown(
                f"- **{row['Database']}**: {row['Link']}",
                unsafe_allow_html=True,
            )
