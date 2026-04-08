"""Chemical Space page - dimensionality reduction and interactive scatter plot."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from chemfuse.web._utils import _find_smiles_column, _params_changed


def render() -> None:
    """Render the Chemical Space page."""
    st.header("Chemical Space Visualization")
    st.caption("Visualize compound distributions using dimensionality reduction.")

    # Data source selection
    source_choice = st.radio(
        "Data source",
        ["Session results", "Upload CSV"],
        horizontal=True,
        key="chemspace_source",
    )

    df: pd.DataFrame | None = None
    _upload_name: str | None = None

    if source_choice == "Session results":
        compounds = st.session_state.get("session_compounds", [])
        if compounds:
            df = pd.DataFrame(compounds)
            st.success(f"Using {len(df)} compounds from session.")
        else:
            st.info("No compounds in session. Run a search first or upload a CSV.")

    else:
        uploaded = st.file_uploader(
            "Upload CSV with SMILES column",
            type=["csv"],
            key="chemspace_upload",
        )
        if uploaded is not None:
            _upload_name = uploaded.name
            try:
                df = pd.read_csv(uploaded, nrows=10_000)
                st.success(f"Loaded {len(df)} rows.")
            except Exception as exc:
                st.error(f"Failed to read CSV: {exc}")

    if df is None or df.empty:
        return

    # Find SMILES column
    smiles_col = _find_smiles_column(df.columns)
    if smiles_col is None:
        st.error("No 'smiles' column found.")
        return

    # Analysis options
    col_method, col_fp = st.columns(2)
    with col_method:
        method = st.selectbox(
            "Dimensionality reduction",
            ["PCA", "UMAP", "t-SNE"],
            key="chemspace_method",
        )
    with col_fp:
        fp_type = st.selectbox(
            "Fingerprint type",
            ["morgan", "maccs", "rdkit"],
            key="chemspace_fp_type",
        )

    # Color-by option
    numeric_cols = [c for c in df.columns if c != smiles_col and pd.to_numeric(df[c], errors="coerce").notna().any()]
    color_options = ["cluster"] + numeric_cols
    color_by = st.selectbox(
        "Color by",
        color_options,
        key="chemspace_color_by",
    )

    # Detect parameter changes and clear stale plot data when they occur.
    _current_params = {
        "source": source_choice,
        "upload_name": _upload_name,
        "method": method,
        "fp_type": fp_type,
        "color_by": color_by,
    }
    if _params_changed(st.session_state, "_chemspace_params_hash", _current_params):
        st.session_state.pop("chemspace_plot_data", None)

    if st.button("Compute Chemical Space", type="primary", key="chemspace_compute_btn"):
        _compute_and_plot(df, smiles_col, method, fp_type, color_by)

    # Show cached plot
    if "chemspace_plot_data" in st.session_state:
        _render_plot(st.session_state["chemspace_plot_data"], color_by)



def _compute_and_plot(
    df: pd.DataFrame,
    smiles_col: str,
    method: str,
    fp_type: str,
    color_by: str,
) -> None:
    """Compute dimensionality reduction and prepare plot data."""
    smiles_list = df[smiles_col].dropna().tolist()
    if len(smiles_list) < 3:
        st.warning("At least 3 valid SMILES required for visualization.")
        return

    with st.spinner(f"Computing {method} projection..."):
        try:
            from chemfuse.analyze.chemspace import reduce_dimensions
            coords = reduce_dimensions(smiles_list, method=method.lower(), fp_type=fp_type)
        except Exception as exc:
            st.error(f"Chemical space computation failed: {exc}")
            return

    if coords is None or len(coords) == 0:
        st.error("Could not compute chemical space coordinates.")
        return

    import numpy as np
    coords_array = np.array(coords)
    plot_df = df.iloc[: len(coords_array)].copy()
    plot_df["x"] = coords_array[:, 0]
    plot_df["y"] = coords_array[:, 1]

    # Add cluster labels if needed
    if color_by == "cluster" and "cluster" not in plot_df.columns:
        try:
            from chemfuse.analyze.clustering import cluster_compounds
            labels = cluster_compounds(smiles_list)
            plot_df["cluster"] = list(labels)[: len(plot_df)]
        except Exception:
            plot_df["cluster"] = 0

    st.session_state["chemspace_plot_data"] = plot_df
    st.session_state["chemspace_method"] = method
    _render_plot(plot_df, color_by)


def _render_plot(df: pd.DataFrame, color_by: str) -> None:
    """Render interactive Plotly scatter plot."""
    try:
        import plotly.express as px
    except ImportError:
        st.error("Plotly required: pip install plotly")
        return

    if "x" not in df.columns or "y" not in df.columns:
        st.warning("No plot data available.")
        return

    method = st.session_state.get("chemspace_method", "PCA")

    # Hover data
    hover_cols: list[str] = []
    for col in ["name", "cid", "smiles", "molecular_weight", "xlogp"]:
        if col in df.columns:
            hover_cols.append(col)

    color_col = color_by if color_by in df.columns else None

    try:
        fig = px.scatter(
            df,
            x="x",
            y="y",
            color=color_col,
            hover_data=hover_cols if hover_cols else None,
            title=f"Chemical Space ({method})",
            labels={"x": f"{method} 1", "y": f"{method} 2"},
        )
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)
    except Exception as exc:
        st.error(f"Plot rendering failed: {exc}")
