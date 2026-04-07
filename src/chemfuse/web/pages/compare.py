"""Compare page - side-by-side compound comparison."""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

_PROPERTIES_TO_COMPARE = [
    ("molecular_weight", "MW (Da)"),
    ("xlogp", "XLogP"),
    ("tpsa", "TPSA (A²)"),
    ("hbd_count", "HBD"),
    ("hba_count", "HBA"),
    ("rotatable_bonds", "Rotatable Bonds"),
    ("heavy_atom_count", "Heavy Atoms"),
    ("complexity", "Complexity"),
]

_DRUGLIKENESS_FILTERS = ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]


def _stable_id(compound: dict[str, Any]) -> str:
    """Return a stable identifier for a compound dict.

    Priority: CID (int) -> SMILES -> name.  The returned string is suitable
    as a dictionary key and does not change when the surrounding list grows or
    shrinks.
    """
    cid = compound.get("cid")
    if cid is not None:
        return f"cid:{cid}"
    smiles = compound.get("smiles")
    if smiles:
        return f"smiles:{smiles}"
    name = compound.get("name")
    if name:
        return f"name:{name}"
    return f"idx:{id(compound)}"


def _display_label(compound: dict[str, Any]) -> str:
    """Return a human-readable label for a compound dict."""
    cid = compound.get("cid")
    name = compound.get("name")
    smiles = compound.get("smiles", "")
    if cid and name:
        return f"CID: {cid} - {name}"
    if cid:
        return f"CID: {cid}"
    if name:
        return name
    return smiles[:30] if smiles else "Unknown"


def render() -> None:
    """Render the Compare page."""
    st.header("Compound Comparison")
    st.caption("Select 2–4 compounds from your session for side-by-side comparison.")

    session_compounds = st.session_state.get("session_compounds", [])
    if len(session_compounds) < 2:
        st.info(
            "Not enough compounds in session. "
            "Run a search to populate the compound pool (need at least 2)."
        )
        return

    # Build stable identifiers for each compound so that adding or removing
    # other compounds from the session does not shift the selection.
    # Priority: CID -> SMILES -> name.  Duplicates are disambiguated with a
    # numeric suffix so the lookup dict remains injective.
    id_to_compound: dict[str, dict[str, Any]] = {}
    id_to_label: dict[str, str] = {}
    for c in session_compounds:
        stable_id = _stable_id(c)
        # Deduplicate: if the same key already exists append a counter
        base = stable_id
        suffix = 1
        while stable_id in id_to_compound:
            stable_id = f"{base}#{suffix}"
            suffix += 1
        id_to_compound[stable_id] = c
        id_to_label[stable_id] = _display_label(c)

    selected_ids = st.multiselect(
        "Select compounds to compare (2–4)",
        options=list(id_to_compound.keys()),
        format_func=lambda sid: id_to_label[sid],
        key="compare_selection",
        max_selections=4,
    )

    if len(selected_ids) < 2:
        st.info("Select at least 2 compounds.")
        return

    selected = [id_to_compound[sid] for sid in selected_ids]
    names = [id_to_label[sid] for sid in selected_ids]

    # Render comparison sections
    _render_property_comparison(selected, names)
    _render_druglikeness_comparison(selected, names)
    _render_admet_comparison(selected, names)


def _render_property_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render side-by-side property comparison table and radar chart."""
    st.subheader("Property Comparison")

    rows: list[dict[str, Any]] = []
    chart_data: dict[str, list[float | None]] = {}

    for prop_key, prop_label in _PROPERTIES_TO_COMPARE:
        row: dict[str, Any] = {"Property": prop_label}
        chart_data[prop_label] = []
        for name, compound in zip(names, compounds, strict=False):
            val = _get_property(compound, prop_key)
            row[name] = f"{val:.2f}" if isinstance(val, float) else (str(val) if val is not None else "—")
            chart_data[prop_label].append(val if isinstance(val, (int, float)) else None)
        rows.append(row)

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

    # Radar chart overlay
    _render_radar_overlay(chart_data, names)


def _render_radar_overlay(
    chart_data: dict[str, list[float | None]],
    names: list[str],
) -> None:
    """Render radar chart overlay for all selected compounds."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        st.caption("Install plotly for radar chart visualization.")
        return

    # Filter to properties with at least one non-None value (allow negative values like XLogP)
    valid_props = {
        k: v for k, v in chart_data.items()
        if any(x is not None for x in v)
    }
    if len(valid_props) < 3:
        return

    categories = list(valid_props.keys())
    fig = go.Figure()

    for i, name in enumerate(names):
        values = [valid_props[prop][i] or 0 for prop in categories]
        values += values[:1]  # Close the polygon
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=categories + [categories[0]],
            name=name,
            fill="toself",
            opacity=0.4,
        ))

    fig.update_layout(
        polar={"radialaxis": {"visible": True}},
        title="Property Radar Chart",
        height=450,
    )
    st.plotly_chart(fig, use_container_width=True)


def _render_druglikeness_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render drug-likeness pass/fail comparison table."""
    st.subheader("Drug-likeness Comparison")

    rows: list[dict[str, Any]] = []
    for filter_name in _DRUGLIKENESS_FILTERS:
        row: dict[str, Any] = {"Filter": filter_name}
        for name, compound in zip(names, compounds, strict=False):
            passed = _get_druglikeness_result(compound, filter_name)
            if passed is None:
                row[name] = "—"
            elif passed:
                row[name] = "PASS"
            else:
                row[name] = "FAIL"
        rows.append(row)

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)


def _render_admet_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render ADMET comparison if available for selected compounds."""
    # Check if any compound has ADMET data
    admet_keys: set[str] = set()
    for compound in compounds:
        dl = compound.get("druglikeness")
        if isinstance(dl, dict):
            admet_keys.update(k for k in dl if k.startswith("admet_"))
        # Also check dedicated ADMET fields
        for admet_field in ("admet", "_admet_profile"):
            admet_data = compound.get(admet_field)
            if isinstance(admet_data, dict):
                admet_keys.update(admet_data.keys())

    if not admet_keys:
        with st.expander("ADMET Comparison"):
            st.caption("Run ADMET predictions on compounds first (via Batch Screen or Profile page).")
        return

    with st.expander("ADMET Comparison"):
        rows = []
        for key in sorted(admet_keys):
            row: dict[str, Any] = {"ADMET Property": key}
            for name, compound in zip(names, compounds, strict=False):
                val = None
                dl = compound.get("druglikeness", {})
                if isinstance(dl, dict):
                    val = dl.get(key)
                if val is None:
                    for admet_field in ("admet", "_admet_profile"):
                        admet_data = compound.get(admet_field)
                        if isinstance(admet_data, dict):
                            val = admet_data.get(key)
                            if val is not None:
                                break
                row[name] = f"{val:.3f}" if isinstance(val, float) else str(val) if val is not None else "—"
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def _get_property(compound: dict[str, Any], key: str) -> Any:
    """Extract a property value from a compound dict, checking nested properties."""
    val = compound.get(key)
    if val is None:
        props = compound.get("properties", {})
        if isinstance(props, dict):
            val = props.get(key)
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return val


def _get_druglikeness_result(compound: dict[str, Any], filter_name: str) -> bool | None:
    """Extract pass/fail result for a drug-likeness filter."""
    dl = compound.get("druglikeness")
    if dl is None:
        return None

    if isinstance(dl, dict):
        filters = dl.get("filters", {})
        if isinstance(filters, dict):
            result = filters.get(filter_name)
            if isinstance(result, dict):
                return result.get("pass")
            if isinstance(result, bool):
                return result

    if hasattr(dl, "filters"):
        filters = dl.filters
        if isinstance(filters, dict):
            result = filters.get(filter_name)
            if isinstance(result, dict):
                return result.get("pass")

    return None
