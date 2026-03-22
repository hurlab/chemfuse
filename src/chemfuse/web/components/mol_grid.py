"""Molecule grid component - displays a grid of compound cards with 2D structure images.

Uses RDKit SVG generation when available, falls back to PubChem PNG URLs.
"""

from __future__ import annotations

import logging
import re
from typing import Any

import streamlit as st

logger = logging.getLogger(__name__)

# Try importing RDKit for SVG generation
_RDKIT_AVAILABLE = False
try:
    from rdkit import Chem  # noqa: F401

    _RDKIT_AVAILABLE = True
except ImportError:
    pass


def _sanitize_svg(svg: str) -> str:
    """Remove script tags and event handlers from SVG to prevent XSS.

    Args:
        svg: Raw SVG string from RDKit.

    Returns:
        Sanitized SVG string safe for embedding with unsafe_allow_html.
    """
    # Remove script elements entirely
    svg = re.sub(r'<script[^>]*>.*?</script>', '', svg, flags=re.DOTALL | re.IGNORECASE)
    # Remove on* event handler attributes (double-quoted values)
    svg = re.sub(r'\s+on\w+\s*=\s*"[^"]*"', '', svg, flags=re.IGNORECASE)
    # Remove on* event handler attributes (single-quoted values)
    svg = re.sub(r"\s+on\w+\s*=\s*'[^']*'", '', svg, flags=re.IGNORECASE)
    return svg


def is_rdkit_available() -> bool:
    """Return True if RDKit is installed and available."""
    return _RDKIT_AVAILABLE


def smiles_to_svg(smiles: str, width: int = 300, height: int = 200) -> str | None:
    """Generate an SVG depiction of a molecule from SMILES using RDKit.

    Args:
        smiles: SMILES string.
        width: Image width in pixels.
        height: Image height in pixels.

    Returns:
        SVG string or None if generation fails.
    """
    if not _RDKIT_AVAILABLE or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        from rdkit.Chem.Draw import rdMolDraw2D

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as exc:
        logger.debug("RDKit SVG generation failed: %s", exc)
        return None


def get_pubchem_png_url(cid: int | None) -> str | None:
    """Get a PubChem 2D structure PNG URL for a compound CID.

    Args:
        cid: PubChem compound ID.

    Returns:
        URL string or None if no CID.
    """
    if cid is None:
        return None
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"


def render_mol_card(
    compound: dict[str, Any],
    width: int = 200,
    height: int = 150,
    on_click_session_key: str | None = None,
) -> None:
    """Render a single molecule card with 2D structure image and basic info.

    Args:
        compound: Compound data dictionary.
        width: Card image width.
        height: Card image height.
        on_click_session_key: Session state key to set when the compound is selected.
    """
    smiles = compound.get("smiles", "")
    cid = compound.get("cid")
    name = compound.get("name") or compound.get("smiles", "Unknown")[:20]

    # Try RDKit SVG first, then PubChem PNG fallback
    svg_data = None
    if smiles and _RDKIT_AVAILABLE:
        svg_data = smiles_to_svg(smiles, width=width, height=height)

    if svg_data:
        safe_svg = _sanitize_svg(svg_data)
        st.markdown(
            f'<div class="mol-card">{safe_svg}</div>',
            unsafe_allow_html=True,
        )
    elif cid:
        img_url = get_pubchem_png_url(cid)
        st.image(img_url, use_container_width=True)
    else:
        st.markdown(
            f'<div class="mol-card" style="padding:40px;color:#888">No structure available<br>'
            f'<code>{smiles[:30] if smiles else "—"}</code></div>',
            unsafe_allow_html=True,
        )

    st.caption(name)

    # Select button
    if on_click_session_key and st.button(
        "View Profile",
        key=f"mol_select_{name}_{cid or hash(smiles or '') or 'no_id'}",
    ):
        st.session_state[on_click_session_key] = compound
        st.session_state["current_page"] = "Profile"
        st.rerun()


def render_mol_grid(
    compounds: list[dict[str, Any]],
    columns: int = 4,
    max_compounds: int = 50,
    on_click_session_key: str | None = None,
) -> None:
    """Render a responsive grid of molecule cards.

    Args:
        compounds: List of compound data dictionaries.
        columns: Number of grid columns.
        max_compounds: Maximum number of compounds to display.
        on_click_session_key: Session state key for compound selection.
    """
    if not compounds:
        st.info("No compounds to display.")
        return

    display = compounds[:max_compounds]
    if len(compounds) > max_compounds:
        st.caption(f"Showing {max_compounds} of {len(compounds)} compounds.")

    cols = st.columns(columns)
    for i, compound in enumerate(display):
        with cols[i % columns]:
            render_mol_card(
                compound,
                on_click_session_key=on_click_session_key,
            )
