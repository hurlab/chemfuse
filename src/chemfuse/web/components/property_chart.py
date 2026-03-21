"""Property chart component - Plotly radar and bar charts for compound properties."""

from __future__ import annotations

import logging
from typing import Any

import streamlit as st

logger = logging.getLogger(__name__)


def render_radar_chart(
    properties: dict[str, float],
    title: str = "Property Profile",
    max_values: dict[str, float] | None = None,
) -> None:
    """Render a Plotly radar chart for compound properties.

    Args:
        properties: Dictionary mapping property name to numeric value.
        title: Chart title.
        max_values: Optional dictionary of maximum reference values for normalization.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        st.caption("Install plotly for property charts: pip install plotly")
        return

    if len(properties) < 3:
        st.caption("At least 3 properties required for radar chart.")
        return

    categories = list(properties.keys())
    values = list(properties.values())

    # Optional normalization against reference maxima
    if max_values:
        normalized = [
            (v / max_values.get(k, v)) if max_values.get(k, 0) > 0 else 0.0
            for k, v in properties.items()
        ]
    else:
        normalized = values

    # Close the polygon
    categories_closed = categories + [categories[0]]
    values_closed = normalized + [normalized[0]]

    fig = go.Figure(
        data=go.Scatterpolar(
            r=values_closed,
            theta=categories_closed,
            fill="toself",
            name=title,
        )
    )
    fig.update_layout(
        polar={"radialaxis": {"visible": True, "showticklabels": True}},
        title=title,
        height=350,
        showlegend=False,
    )
    st.plotly_chart(fig, use_container_width=True)


def render_bar_chart(
    properties: dict[str, float],
    title: str = "Property Comparison",
    color: str = "#1f77b4",
) -> None:
    """Render a Plotly bar chart for compound properties.

    Args:
        properties: Dictionary mapping property name to numeric value.
        title: Chart title.
        color: Bar color (hex or named color).
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        st.caption("Install plotly for property charts: pip install plotly")
        return

    if not properties:
        return

    fig = go.Figure(
        data=go.Bar(
            x=list(properties.keys()),
            y=list(properties.values()),
            marker_color=color,
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Property",
        yaxis_title="Value",
        height=350,
    )
    st.plotly_chart(fig, use_container_width=True)


def render_multi_radar_chart(
    compounds_properties: list[dict[str, float]],
    names: list[str],
    title: str = "Property Comparison",
) -> None:
    """Render a multi-trace radar chart comparing several compounds.

    Args:
        compounds_properties: List of property dictionaries, one per compound.
        names: Compound names (one per property dict).
        title: Chart title.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        st.caption("Install plotly for property charts: pip install plotly")
        return

    if not compounds_properties or len(compounds_properties) != len(names):
        return

    # Use the union of all property keys
    all_keys = list(
        {k for props in compounds_properties for k in props}
    )
    if len(all_keys) < 3:
        return

    fig = go.Figure()
    for name, props in zip(names, compounds_properties, strict=False):
        values = [props.get(k, 0.0) for k in all_keys]
        # Close the polygon
        theta = all_keys + [all_keys[0]]
        r = values + [values[0]]
        fig.add_trace(go.Scatterpolar(r=r, theta=theta, fill="toself", name=name, opacity=0.5))

    fig.update_layout(
        polar={"radialaxis": {"visible": True}},
        title=title,
        height=450,
    )
    st.plotly_chart(fig, use_container_width=True)


def extract_chart_properties(compound: dict[str, Any]) -> dict[str, float]:
    """Extract plottable numeric properties from a compound dictionary.

    Args:
        compound: Compound data dictionary.

    Returns:
        Dictionary of property name to float value for charting.
    """
    prop_keys = [
        ("molecular_weight", "MW"),
        ("xlogp", "XLogP"),
        ("tpsa", "TPSA"),
        ("hbd_count", "HBD"),
        ("hba_count", "HBA"),
        ("rotatable_bonds", "RotBonds"),
    ]
    result: dict[str, float] = {}
    for key, label in prop_keys:
        val = compound.get(key)
        if val is None:
            props = compound.get("properties", {})
            if isinstance(props, dict):
                val = props.get(key)
        if val is not None:
            try:
                result[label] = float(val)
            except (ValueError, TypeError):
                pass
    return result
