"""Batch Screen page - CSV upload and pipeline execution."""

from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st


def render() -> None:
    """Render the Batch Screen page."""
    st.header("Batch Screening")
    st.caption("Upload a CSV file with a SMILES column to screen a compound library.")

    # File upload
    uploaded_file = st.file_uploader(
        "Upload CSV",
        type=["csv"],
        key="batch_upload",
        help="CSV must contain a column named 'smiles' (case-insensitive).",
    )

    if uploaded_file is None:
        _show_example_format()
        return

    # Parse uploaded file
    try:
        df_input = pd.read_csv(uploaded_file)
    except Exception as exc:
        st.error(f"Failed to read CSV: {exc}")
        return

    smiles_col = _find_smiles_column(df_input)
    if smiles_col is None:
        st.error("No 'smiles' column found in the uploaded CSV.")
        return

    st.success(f"Loaded {len(df_input)} rows. SMILES column: '{smiles_col}'")

    # Pipeline configuration
    st.subheader("Pipeline Configuration")
    col_src, col_opts = st.columns(2)

    with col_src:
        st.markdown("**Sources**")
        use_pubchem = st.checkbox("PubChem enrichment", value=False, key="batch_pubchem")
        use_chembl = st.checkbox("ChEMBL enrichment", value=False, key="batch_chembl")

    with col_opts:
        st.markdown("**Analysis Options**")
        run_admet = st.toggle("ADMET predictions", value=False, key="batch_admet")
        run_druglikeness = st.toggle("Drug-likeness filters", value=True, key="batch_dl")
        run_clustering = st.toggle("Clustering", value=False, key="batch_cluster")

    # Run button
    if st.button("Run Screening Pipeline", type="primary", key="batch_run_btn"):
        sources = []
        if use_pubchem:
            sources.append("pubchem")
        if use_chembl:
            sources.append("chembl")
        _run_pipeline(
            df_input,
            smiles_col,
            sources,
            run_admet=run_admet,
            run_druglikeness=run_druglikeness,
            run_clustering=run_clustering,
        )

    # Display results
    results_df = st.session_state.get("batch_results")
    if results_df is not None:
        _display_results(results_df)


def _show_example_format() -> None:
    """Show an example of the expected CSV format."""
    with st.expander("Expected CSV format"):
        example = pd.DataFrame({
            "smiles": ["CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
            "name": ["aspirin", "caffeine"],
        })
        st.dataframe(example, use_container_width=True, hide_index=True)


def _find_smiles_column(df: pd.DataFrame) -> str | None:
    """Find the SMILES column in a dataframe (case-insensitive)."""
    for col in df.columns:
        if col.lower() == "smiles":
            return col
    return None


def _run_pipeline(
    df: pd.DataFrame,
    smiles_col: str,
    sources: list[str],
    *,
    run_admet: bool,
    run_druglikeness: bool,
    run_clustering: bool,
) -> None:
    """Execute the screening pipeline with progress reporting."""
    smiles_list = df[smiles_col].dropna().tolist()
    if not smiles_list:
        st.error("No valid SMILES found in the uploaded file.")
        return

    total_steps = 1 + (1 if run_druglikeness else 0) + (1 if run_admet else 0) + (1 if run_clustering else 0)
    progress = st.progress(0, text="Initializing pipeline...")
    step = 0

    # Step 1: Build base dataframe
    results = df.copy()
    step += 1
    progress.progress(step / total_steps, text="Loading compounds...")

    # Step 2: Drug-likeness
    if run_druglikeness:
        progress.progress(step / total_steps, text="Computing drug-likeness...")
        dl_results = _compute_druglikeness(smiles_list)
        for key, values in dl_results.items():
            results[key] = values
        step += 1
        progress.progress(step / total_steps, text="Drug-likeness complete.")

    # Step 3: ADMET predictions
    if run_admet:
        progress.progress(step / total_steps, text="Running ADMET predictions...")
        admet_results = _compute_admet(smiles_list)
        for key, values in admet_results.items():
            results[f"admet_{key}"] = values
        step += 1
        progress.progress(step / total_steps, text="ADMET complete.")

    # Step 4: Clustering
    if run_clustering:
        progress.progress(step / total_steps, text="Clustering compounds...")
        cluster_labels = _compute_clusters(smiles_list)
        results["cluster"] = cluster_labels
        step += 1
        progress.progress(step / total_steps, text="Clustering complete.")

    progress.progress(1.0, text="Pipeline complete!")

    st.session_state["batch_results"] = results
    st.success(f"Screening complete: {len(results)} compounds processed.")


def _compute_druglikeness(smiles_list: list[str]) -> dict[str, list[Any]]:
    """Compute drug-likeness for a list of SMILES."""
    lipinski_pass = []
    mw_values = []

    for smi in smiles_list:
        try:
            from chemfuse.compute.druglikeness import check_drug_likeness
            dl = check_drug_likeness({}, smiles=smi)
            if hasattr(dl, "filters") and dl.filters:
                lip = dl.filters.get("Lipinski", {})
                passed = lip.get("pass", False) if isinstance(lip, dict) else False
                lipinski_pass.append(passed)
            else:
                lipinski_pass.append(None)
        except Exception:
            lipinski_pass.append(None)
        mw_values.append(None)

    return {"lipinski_pass": lipinski_pass}


def _compute_admet(smiles_list: list[str]) -> dict[str, list[Any]]:
    """Compute ADMET predictions for a list of SMILES."""
    try:
        from chemfuse.compute.admet import predict_admet
        preds = predict_admet(smiles_list)
        if preds and isinstance(preds, list) and preds:
            keys = list(preds[0].keys()) if isinstance(preds[0], dict) else []
            return {k: [p.get(k) for p in preds] for k in keys}
    except Exception:
        pass
    return {}


def _compute_clusters(smiles_list: list[str]) -> list[int | None]:
    """Compute cluster labels for a list of SMILES."""
    try:
        from chemfuse.analyze.clustering import cluster_compounds
        labels = cluster_compounds(smiles_list)
        return list(labels)
    except Exception:
        return [None] * len(smiles_list)


def _display_results(df: pd.DataFrame) -> None:
    """Display batch screening results with download option."""
    st.subheader(f"Screening Results ({len(df)} compounds)")

    # Column filters
    with st.expander("Filter Results"):
        from chemfuse.web.components.filters import render_property_filters
        filter_params = render_property_filters()
        if filter_params:
            df = _apply_filters(df, filter_params)

    st.dataframe(df, use_container_width=True, hide_index=True)

    # Download button
    col_csv, col_xlsx = st.columns(2)
    with col_csv:
        buf = io.BytesIO()
        df.to_csv(buf, index=False)
        buf.seek(0)
        st.download_button(
            "Download CSV",
            data=buf,
            file_name="chemfuse_screen_results.csv",
            mime="text/csv",
            key="batch_download_csv",
        )
    with col_xlsx:
        try:
            buf_xl = io.BytesIO()
            df.to_excel(buf_xl, index=False)
            buf_xl.seek(0)
            st.download_button(
                "Download Excel",
                data=buf_xl,
                file_name="chemfuse_screen_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="batch_download_xlsx",
            )
        except ImportError:
            st.caption("Excel export requires: pip install openpyxl")


def _apply_filters(df: pd.DataFrame, filters: dict[str, tuple[float, float]]) -> pd.DataFrame:
    """Apply numeric range filters to a dataframe."""
    mask = pd.Series([True] * len(df), index=df.index)
    for col, (min_val, max_val) in filters.items():
        if col in df.columns:
            numeric = pd.to_numeric(df[col], errors="coerce")
            mask &= (numeric >= min_val) & (numeric <= max_val)
    return df[mask].reset_index(drop=True)
