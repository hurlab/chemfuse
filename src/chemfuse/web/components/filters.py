"""Sidebar filter widgets for compound property ranges and drug-likeness."""

from __future__ import annotations

import streamlit as st

# Default property ranges for common physicochemical properties
_PROPERTY_DEFAULTS: dict[str, tuple[float, float, float, float]] = {
    # key: (absolute_min, absolute_max, default_min, default_max)
    "molecular_weight": (0.0, 2000.0, 0.0, 1000.0),
    "xlogp": (-20.0, 20.0, -5.0, 10.0),
    "tpsa": (0.0, 500.0, 0.0, 200.0),
    "hbd_count": (0.0, 30.0, 0.0, 15.0),
    "hba_count": (0.0, 30.0, 0.0, 15.0),
    "rotatable_bonds": (0.0, 50.0, 0.0, 20.0),
}

_PROPERTY_LABELS: dict[str, str] = {
    "molecular_weight": "MW (Da)",
    "xlogp": "XLogP",
    "tpsa": "TPSA (A²)",
    "hbd_count": "HBD",
    "hba_count": "HBA",
    "rotatable_bonds": "Rotatable Bonds",
}

_DRUGLIKENESS_FILTER_NAMES = ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]


def render_property_filters(
    show_mw: bool = True,
    show_logp: bool = True,
    show_tpsa: bool = True,
    show_hbd: bool = False,
    show_hba: bool = False,
    show_rotbonds: bool = False,
) -> dict[str, tuple[float, float]]:
    """Render sidebar range filter widgets and return active filter ranges.

    Args:
        show_mw: Include molecular weight filter.
        show_logp: Include XLogP filter.
        show_tpsa: Include TPSA filter.
        show_hbd: Include HBD filter.
        show_hba: Include HBA filter.
        show_rotbonds: Include rotatable bonds filter.

    Returns:
        Dictionary mapping property key to (min_val, max_val) tuple.
    """
    active_filters: dict[str, tuple[float, float]] = {}

    filter_flags = {
        "molecular_weight": show_mw,
        "xlogp": show_logp,
        "tpsa": show_tpsa,
        "hbd_count": show_hbd,
        "hba_count": show_hba,
        "rotatable_bonds": show_rotbonds,
    }

    for prop_key, enabled in filter_flags.items():
        if not enabled:
            continue
        abs_min, abs_max, def_min, def_max = _PROPERTY_DEFAULTS[prop_key]
        label = _PROPERTY_LABELS[prop_key]
        step = 1.0 if prop_key.endswith("_count") or prop_key == "rotatable_bonds" else 0.1

        values = st.slider(
            label,
            min_value=float(abs_min),
            max_value=float(abs_max),
            value=(float(def_min), float(def_max)),
            step=step,
            key=f"filter_{prop_key}",
        )
        active_filters[prop_key] = values

    return active_filters


def render_druglikeness_checkboxes() -> list[str]:
    """Render drug-likeness filter checkboxes and return selected filter names.

    Returns:
        List of selected drug-likeness filter names.
    """
    st.markdown("**Drug-likeness Filters**")
    selected = []
    for name in _DRUGLIKENESS_FILTER_NAMES:
        if st.checkbox(name, key=f"dl_filter_{name}"):
            selected.append(name)
    return selected


def render_full_sidebar_filters() -> dict[str, object]:
    """Render a complete sidebar filter panel including property ranges and drug-likeness.

    Returns:
        Dictionary with keys 'property_ranges' and 'druglikeness_filters'.
    """
    st.markdown("**Property Filters**")
    property_ranges = render_property_filters(
        show_mw=True,
        show_logp=True,
        show_tpsa=True,
        show_hbd=True,
        show_hba=True,
        show_rotbonds=True,
    )

    st.divider()
    druglikeness_filters = render_druglikeness_checkboxes()

    return {
        "property_ranges": property_ranges,
        "druglikeness_filters": druglikeness_filters,
    }
