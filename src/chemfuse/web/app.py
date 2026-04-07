"""ChemFuse Web UI - Main Streamlit application entry point.

Multi-page Streamlit app with sidebar navigation for ChemFuse.
Run with: streamlit run -m chemfuse.web.app
or via: chemfuse web
"""

from __future__ import annotations

import streamlit as st

from chemfuse._version import __version__

# CSS for molecule card styling
_MOL_CARD_CSS = """
<style>
.mol-card {
    border: 1px solid #e0e0e0;
    border-radius: 8px;
    padding: 8px;
    text-align: center;
    background: #fafafa;
    margin-bottom: 8px;
}
.mol-card:hover {
    border-color: #1f77b4;
    background: #f0f7ff;
}
.pass-badge {
    color: #2e7d32;
    font-weight: bold;
}
.fail-badge {
    color: #c62828;
    font-weight: bold;
}
</style>
"""

_PAGES: dict[str, str] = {
    "Search": "Search compounds across PubChem and ChEMBL",
    "Profile": "Deep view of a single compound",
    "Batch Screen": "Screen compound libraries from CSV",
    "Chemical Space": "Visualize chemical space with dimensionality reduction",
    "Compare": "Side-by-side comparison of 2-4 compounds",
    "Cross-Reference": "Map identifiers across databases",
}


def _init_session_state() -> None:
    """Initialize session state keys on first load."""
    defaults: dict[str, object] = {
        "current_page": "Search",
        "search_results": [],
        "selected_compound": None,
        "session_compounds": [],
        "batch_results": None,
        "xref_result": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def _render_sidebar() -> str:
    """Render sidebar navigation and return the selected page name."""
    with st.sidebar:
        st.title("ChemFuse")
        st.caption(f"v{__version__} — Multi-database cheminformatics suite")
        st.divider()

        selected = st.radio(
            "Navigation",
            list(_PAGES.keys()),
            index=list(_PAGES.keys()).index(st.session_state.current_page),
            key="nav_radio",
        )
        st.session_state.current_page = selected
        st.caption(_PAGES[selected])
        st.divider()
        st.markdown("**Data Sources**")
        st.caption("PubChem • ChEMBL • BindingDB\nSureChEMBL • Open Targets • UniChem")
        st.divider()
        if st.button("Clear Session", key="clear_session"):
            st.session_state.clear()
            st.rerun()

    return selected


def _render_page(page_name: str) -> None:
    """Dispatch rendering to the appropriate page module."""
    try:
        if page_name == "Search":
            from chemfuse.web.pages.search import render
            render()
        elif page_name == "Profile":
            from chemfuse.web.pages.profile import render
            render()
        elif page_name == "Batch Screen":
            from chemfuse.web.pages.screen import render
            render()
        elif page_name == "Chemical Space":
            from chemfuse.web.pages.chemspace import render
            render()
        elif page_name == "Compare":
            from chemfuse.web.pages.compare import render
            render()
        elif page_name == "Cross-Reference":
            from chemfuse.web.pages.xref import render
            render()
        else:
            st.error(f"Unknown page: {page_name}")
    except Exception as exc:
        st.error(f"Error loading page '{page_name}': {exc}")


def main() -> None:
    """Main Streamlit application entry point."""
    st.set_page_config(
        page_title="ChemFuse",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.markdown(_MOL_CARD_CSS, unsafe_allow_html=True)

    _init_session_state()
    selected_page = _render_sidebar()
    _render_page(selected_page)


if __name__ == "__main__":
    main()
