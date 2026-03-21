"""Cross-Reference page - identifier mapping across databases."""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

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
            result = _lookup_xref(identifier.strip())
        if result:
            st.session_state["xref_result"] = result
        else:
            st.warning("No cross-references found for this identifier.")
            st.session_state["xref_result"] = None

    # Display result
    xref_data = st.session_state.get("xref_result")
    if xref_data:
        _display_xref_table(identifier.strip(), xref_data)


@st.cache_data(ttl=600, show_spinner=False)
def _lookup_xref(identifier: str) -> dict[str, str] | None:
    """Look up cross-references for an identifier via UniChem."""
    import asyncio

    async def _run() -> dict[str, str] | None:
        try:
            from chemfuse.sources.unichem import UniChemSource
            src = UniChemSource()
            async with src:
                # Determine identifier type
                if identifier.isdigit():
                    # PubChem CID
                    return await src.map_identifiers(identifier, "pubchem")
                elif identifier.upper().startswith("CHEMBL"):
                    return await src.map_identifiers(identifier, "chembl")
                elif len(identifier) == 27 and "-" in identifier:
                    # InChIKey format: XXXXX-XXXXX-X
                    return await src.cross_reference(identifier)
                else:
                    # Try as SMILES or name via PubChem first
                    return await _resolve_via_pubchem(identifier, src)
        except Exception as exc:
            st.error(f"Cross-reference lookup failed: {exc}")
            return None

    loop = asyncio.new_event_loop()
    result = loop.run_until_complete(_run())
    loop.close()
    return result


async def _resolve_via_pubchem(identifier: str, unichem_src: Any) -> dict[str, str] | None:
    """Resolve identifier via PubChem then get cross-references."""
    try:
        from chemfuse.sources.pubchem import PubChemSource
        pubchem_src = PubChemSource()
        async with pubchem_src:
            # Try name search
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
            link = url_template.format(id=db_id)
            link_html = f'<a href="{link}" target="_blank">{db_id}</a>'
        else:
            link_html = db_id
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
