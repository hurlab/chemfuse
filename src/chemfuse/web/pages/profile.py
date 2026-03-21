"""Compound Profile page - deep view of a single compound."""

from __future__ import annotations

from typing import Any

import streamlit as st


def render() -> None:
    """Render the Compound Profile page."""
    st.header("Compound Profile")

    # Get compound from session state
    compound_data = st.session_state.get("selected_compound")

    if compound_data is None:
        st.info(
            "No compound selected. Search for a compound and click on it, "
            "or enter a SMILES/CID below."
        )
        _render_manual_lookup()
        return

    _render_compound_profile(compound_data)


def _render_manual_lookup() -> None:
    """Allow manual compound lookup by SMILES or CID."""
    col_id, col_btn = st.columns([3, 1])
    with col_id:
        identifier = st.text_input(
            "Enter CID or SMILES",
            key="profile_manual_id",
            placeholder="e.g., 2244 or CC(=O)Oc1ccccc1C(=O)O",
        )
    with col_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Look up", key="profile_lookup_btn") and identifier.strip():
            data = _lookup_compound(identifier.strip())
            if data:
                st.session_state["selected_compound"] = data
                st.rerun()
            else:
                st.warning("Compound not found.")


def _lookup_compound(identifier: str) -> dict[str, Any] | None:
    """Fetch compound data by CID or SMILES."""
    import asyncio

    async def _run() -> dict[str, Any] | None:
        try:
            from chemfuse.sources.pubchem import PubChemSource
            src = PubChemSource()
            async with src:
                if identifier.isdigit():
                    results = await src.search(identifier, query_type="cid")
                else:
                    results = await src.search(identifier, query_type="smiles")
            if results:
                c = results[0]
                return c.to_dict() if hasattr(c, "to_dict") else {}
        except Exception as exc:
            st.error(f"Lookup error: {exc}")
        return None

    loop = asyncio.new_event_loop()
    result = loop.run_until_complete(_run())
    loop.close()
    return result


def _render_compound_profile(data: dict[str, Any]) -> None:
    """Render the full compound profile."""
    name = data.get("name") or data.get("smiles", "Unknown")[:30]
    st.subheader(name)

    # Structure image
    col_img, col_ids = st.columns([1, 2])
    with col_img:
        smiles = data.get("smiles", "")
        cid = data.get("cid")
        if cid:
            img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
            st.image(img_url, caption="2D Structure", use_container_width=True)
        elif smiles:
            st.caption(f"SMILES: `{smiles[:60]}`")
        else:
            st.caption("No structure available")

    with col_ids:
        _render_identifiers(data)

    # Properties section
    with st.expander("Physicochemical Properties", expanded=True):
        _render_properties(data)

    # Drug-likeness section
    with st.expander("Drug-likeness Filters", expanded=True):
        _render_druglikeness(data)

    # ADMET predictions
    with st.expander("ADMET Predictions"):
        _render_admet(data)

    # Enrichment sections
    with st.expander("Bioactivities (ChEMBL)"):
        _render_bioactivities(data)

    with st.expander("Binding Data (BindingDB)"):
        _render_binding(data)

    with st.expander("Patents (SureChEMBL)"):
        _render_patents(data)

    with st.expander("Target-Disease Associations (Open Targets)"):
        _render_targets(data)


def _render_identifiers(data: dict[str, Any]) -> None:
    """Render compound identifier section."""
    st.markdown("**Identifiers**")
    rows = []
    for label, key in [
        ("CID", "cid"),
        ("ChEMBL ID", "chembl_id"),
        ("InChIKey", "inchikey"),
        ("Formula", "formula"),
        ("Name", "name"),
    ]:
        val = data.get(key)
        if val:
            rows.append({"Field": label, "Value": str(val)})
    if rows:
        import pandas as pd
        st.dataframe(
            pd.DataFrame(rows),
            use_container_width=True,
            hide_index=True,
        )


def _render_properties(data: dict[str, Any]) -> None:
    """Render physicochemical properties table and radar chart."""
    import pandas as pd

    from chemfuse.web.components.property_chart import render_radar_chart

    prop_map = {
        "molecular_weight": ("Molecular Weight", "Da"),
        "xlogp": ("XLogP", ""),
        "tpsa": ("TPSA", "A²"),
        "hbd_count": ("HBD", ""),
        "hba_count": ("HBA", ""),
        "rotatable_bonds": ("Rotatable Bonds", ""),
        "heavy_atom_count": ("Heavy Atoms", ""),
        "complexity": ("Complexity", ""),
    }
    rows = []
    chart_data: dict[str, float] = {}
    for key, (label, unit) in prop_map.items():
        val = data.get(key)
        if val is None:
            # Try nested properties dict
            props = data.get("properties", {})
            if isinstance(props, dict):
                val = props.get(key)
        if val is not None:
            try:
                float_val = float(val)
                display = f"{float_val:.2f} {unit}".strip()
                rows.append({"Property": label, "Value": display})
                chart_data[label] = float_val
            except (ValueError, TypeError):
                rows.append({"Property": label, "Value": str(val)})

    if rows:
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        if len(chart_data) >= 3:
            render_radar_chart(chart_data, title="Property Profile")


def _render_druglikeness(data: dict[str, Any]) -> None:
    """Render drug-likeness filter pass/fail badges."""
    # Try pre-computed druglikeness data
    dl_data = data.get("druglikeness")

    if dl_data is None:
        # Compute on-the-fly
        smiles = data.get("smiles", "")
        if smiles:
            try:
                from chemfuse.compute.druglikeness import check_drug_likeness
                props = {
                    "molecular_weight": data.get("molecular_weight"),
                    "xlogp": data.get("xlogp"),
                    "tpsa": data.get("tpsa"),
                    "hbd_count": data.get("hbd_count"),
                    "hba_count": data.get("hba_count"),
                    "rotatable_bonds": data.get("rotatable_bonds"),
                }
                dl_data = check_drug_likeness(props, smiles=smiles)
            except Exception:
                pass

    if dl_data is None:
        st.caption("Drug-likeness data not available. SMILES required.")
        return

    # Display filter results
    filters = getattr(dl_data, "filters", None) or (
        dl_data if isinstance(dl_data, dict) else {}
    )
    if hasattr(dl_data, "model_dump"):
        filters = dl_data.model_dump().get("filters", {})

    if not filters:
        st.caption("No drug-likeness filter data available.")
        return

    cols = st.columns(min(len(filters), 5))
    for i, (filter_name, result) in enumerate(filters.items()):
        with cols[i % len(cols)]:
            passed = result.get("pass", False) if isinstance(result, dict) else getattr(result, "pass_", False)
            badge = "PASS" if passed else "FAIL"
            color = "green" if passed else "red"
            st.markdown(
                f"**{filter_name}**  \n"
                f'<span style="color:{color};font-weight:bold">{badge}</span>',
                unsafe_allow_html=True,
            )


def _render_admet(data: dict[str, Any]) -> None:
    """Render ADMET predictions if available."""
    smiles = data.get("smiles", "")
    if not smiles:
        st.caption("SMILES required for ADMET predictions.")
        return

    if st.button("Predict ADMET", key="admet_predict_btn"):
        with st.spinner("Computing ADMET predictions..."):
            try:
                from chemfuse.compute.admet import predict_admet
                preds = predict_admet([smiles])
                if preds:
                    import pandas as pd
                    st.dataframe(
                        pd.DataFrame(preds),
                        use_container_width=True,
                        hide_index=True,
                    )
                else:
                    st.info("No ADMET predictions available.")
            except Exception as exc:
                st.error(f"ADMET prediction failed: {exc}")


def _render_bioactivities(data: dict[str, Any]) -> None:
    """Render ChEMBL bioactivity data with enrich button."""
    bioacts = data.get("bioactivities", [])
    if bioacts:
        import pandas as pd
        st.dataframe(
            pd.DataFrame([b if isinstance(b, dict) else vars(b) for b in bioacts[:50]]),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("No bioactivity data loaded.")

    if st.button("Enrich from ChEMBL", key="enrich_chembl_btn"):
        with st.spinner("Fetching ChEMBL bioactivities..."):
            _enrich_bioactivities(data)


def _enrich_bioactivities(data: dict[str, Any]) -> None:
    """Fetch bioactivities from ChEMBL and update session state."""
    import asyncio

    async def _run() -> list[Any]:
        chembl_id = data.get("chembl_id")
        if not chembl_id:
            st.warning("ChEMBL ID required for bioactivity enrichment.")
            return []
        try:
            from chemfuse.sources.chembl import ChEMBLSource
            src = ChEMBLSource()
            async with src:
                return await src.get_bioactivities(chembl_id)
        except Exception as exc:
            st.error(f"ChEMBL enrichment failed: {exc}")
            return []

    loop = asyncio.new_event_loop()
    bioacts = loop.run_until_complete(_run())
    loop.close()

    if bioacts:
        data["bioactivities"] = bioacts
        st.session_state["selected_compound"] = data
        st.success(f"Loaded {len(bioacts)} bioactivity records.")
        st.rerun()


def _render_binding(data: dict[str, Any]) -> None:
    """Render binding data with enrich button."""
    binding = data.get("binding_data", [])
    if binding:
        import pandas as pd
        st.dataframe(
            pd.DataFrame([b if isinstance(b, dict) else vars(b) for b in binding[:50]]),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("No binding data loaded.")

    if st.button("Enrich from BindingDB", key="enrich_bindingdb_btn"):
        st.info("BindingDB enrichment requires SMILES.")


def _render_patents(data: dict[str, Any]) -> None:
    """Render patent data with enrich button."""
    patents = data.get("patents", [])
    if patents:
        import pandas as pd
        st.dataframe(
            pd.DataFrame([p if isinstance(p, dict) else vars(p) for p in patents[:50]]),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("No patent data loaded.")

    if st.button("Enrich from SureChEMBL", key="enrich_surechembl_btn"):
        st.info("SureChEMBL enrichment requires InChIKey or SMILES.")


def _render_targets(data: dict[str, Any]) -> None:
    """Render target-disease associations with enrich button."""
    targets = data.get("target_associations", [])
    if targets:
        import pandas as pd
        st.dataframe(
            pd.DataFrame([t if isinstance(t, dict) else vars(t) for t in targets[:50]]),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("No target-disease data loaded.")

    if st.button("Enrich from Open Targets", key="enrich_opentargets_btn"):
        st.info("Open Targets enrichment requires ChEMBL ID.")
