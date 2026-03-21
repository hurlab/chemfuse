"""Streamlit mock fixtures for web UI tests.

All tests in this package run WITHOUT a live Streamlit server.
We mock the streamlit module so that rendering functions can be
imported and exercised without triggering actual Streamlit calls.
"""

from __future__ import annotations

import sys
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Minimal session-state stub that supports attribute and key access
# ---------------------------------------------------------------------------

class _SessionState(dict):
    """Dict subclass that also supports attribute-style access."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError:
            return None

    def __setattr__(self, name: str, value: Any) -> None:
        self[name] = value

    def __delattr__(self, name: str) -> None:
        try:
            del self[name]
        except KeyError:
            raise AttributeError(name) from None


@pytest.fixture()
def session_state() -> _SessionState:
    """Return a fresh session state stub."""
    return _SessionState()


# ---------------------------------------------------------------------------
# Streamlit mock builder
# ---------------------------------------------------------------------------

def _make_st_mock(state: _SessionState | None = None) -> MagicMock:
    """Build a comprehensive mock of the streamlit module."""
    st = MagicMock(name="streamlit")

    # session_state
    _state = state if state is not None else _SessionState()
    st.session_state = _state

    # Widgets that return values
    st.text_input.return_value = ""
    st.number_input.return_value = 0
    st.selectbox.return_value = "name"
    st.multiselect.return_value = []
    st.checkbox.return_value = False
    st.toggle.return_value = False
    st.radio.return_value = "Grid"
    st.slider.return_value = (0.0, 1000.0)
    st.file_uploader.return_value = None
    st.button.return_value = False

    # Layout - context managers
    def _make_ctx(*_args: Any, **_kwargs: Any) -> MagicMock:
        ctx = MagicMock()
        ctx.__enter__ = lambda s: s
        ctx.__exit__ = MagicMock(return_value=False)
        return ctx

    # st.columns must return exactly the right number of columns.
    # The argument can be an int (number of columns) or a list of ratios.
    def _columns(spec: Any, *_args: Any, **_kwargs: Any) -> list:
        if isinstance(spec, (list, tuple)):
            n = len(spec)
        else:
            try:
                n = int(spec)
            except (TypeError, ValueError):
                n = 2
        return [_make_ctx() for _ in range(n)]

    st.columns.side_effect = _columns
    st.expander.side_effect = _make_ctx
    st.sidebar = _make_ctx()
    st.container.side_effect = _make_ctx
    st.form.side_effect = _make_ctx

    # Non-returning display calls
    st.write = MagicMock()
    st.markdown = MagicMock()
    st.header = MagicMock()
    st.subheader = MagicMock()
    st.caption = MagicMock()
    st.info = MagicMock()
    st.warning = MagicMock()
    st.error = MagicMock()
    st.success = MagicMock()
    st.divider = MagicMock()
    st.image = MagicMock()
    st.dataframe = MagicMock()
    st.plotly_chart = MagicMock()
    st.progress = MagicMock()

    # spinner must be a MagicMock so tests can assert_called / assert_not_called,
    # but it also needs to work as a context manager (with st.spinner(...): ...)
    _spinner_mock = MagicMock()
    _spinner_mock.side_effect = _make_ctx
    st.spinner = _spinner_mock

    st.download_button = MagicMock()
    st.code = MagicMock()
    st.rerun = MagicMock()
    st.set_page_config = MagicMock()

    # cache_data decorator - passes through the function unchanged
    def _cache_data(*_args: Any, **_kwargs: Any):
        def decorator(fn):
            return fn
        # Handle both @st.cache_data and @st.cache_data(ttl=...)
        if _args and callable(_args[0]):
            return _args[0]
        return decorator

    st.cache_data = _cache_data

    return st


@pytest.fixture()
def mock_st(session_state: _SessionState) -> MagicMock:
    """Provide a pre-configured streamlit mock with a fresh session state."""
    return _make_st_mock(state=session_state)


# ---------------------------------------------------------------------------
# Web module names whose `st` attribute we need to patch
# ---------------------------------------------------------------------------

_WEB_MODULES = [
    "chemfuse.web.pages.search",
    "chemfuse.web.pages.profile",
    "chemfuse.web.pages.screen",
    "chemfuse.web.pages.chemspace",
    "chemfuse.web.pages.compare",
    "chemfuse.web.pages.xref",
    "chemfuse.web.components.mol_grid",
    "chemfuse.web.components.property_chart",
    "chemfuse.web.components.filters",
    "chemfuse.web.app",
]


@pytest.fixture()
def patch_st(mock_st: MagicMock):
    """Patch the streamlit module globally AND in all already-imported web modules."""
    # 1. Patch sys.modules so future imports resolve to the mock.
    # 2. Also patch the `st` attribute in any already-imported web modules so
    #    that modules loaded before this fixture still use the mock.
    patchers: list[Any] = []
    # sys.modules patch
    sys_patch = patch.dict("sys.modules", {"streamlit": mock_st})
    sys_patch.start()
    patchers.append(sys_patch)

    # Per-module `st` attribute patches
    for mod_name in _WEB_MODULES:
        # Only patch modules that are already in sys.modules (i.e. already imported).
        mod = sys.modules.get(mod_name)
        if mod is not None and hasattr(mod, "st"):
            p = patch.object(mod, "st", mock_st)
            p.start()
            patchers.append(p)

    try:
        yield mock_st
    finally:
        for p in reversed(patchers):
            p.stop()


# ---------------------------------------------------------------------------
# Sample compound data fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def aspirin_data() -> dict[str, Any]:
    """Sample aspirin compound data dict."""
    return {
        "cid": 2244,
        "chembl_id": "CHEMBL25",
        "name": "aspirin",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "formula": "C9H8O4",
        "sources": ["pubchem"],
        "molecular_weight": 180.16,
        "xlogp": 1.19,
        "tpsa": 63.6,
        "hbd_count": 1,
        "hba_count": 4,
        "rotatable_bonds": 3,
        "bioactivities": [],
        "binding_data": [],
        "patents": [],
        "target_associations": [],
        "druglikeness": None,
    }


@pytest.fixture()
def caffeine_data() -> dict[str, Any]:
    """Sample caffeine compound data dict."""
    return {
        "cid": 2519,
        "chembl_id": "CHEMBL113",
        "name": "caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "formula": "C8H10N4O2",
        "sources": ["pubchem"],
        "molecular_weight": 194.19,
        "xlogp": -0.07,
        "tpsa": 58.44,
        "hbd_count": 0,
        "hba_count": 6,
        "rotatable_bonds": 0,
        "bioactivities": [],
        "binding_data": [],
        "patents": [],
        "target_associations": [],
        "druglikeness": None,
    }


@pytest.fixture()
def compound_list(aspirin_data: dict, caffeine_data: dict) -> list[dict]:
    """Return a list of sample compound dicts."""
    return [aspirin_data, caffeine_data]
