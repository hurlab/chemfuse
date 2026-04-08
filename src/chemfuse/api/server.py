"""FastAPI REST API for ChemFuse."""

from __future__ import annotations

import logging
from typing import Any

from chemfuse._version import __version__

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pydantic request models (defined at module level for FastAPI compatibility)
# ---------------------------------------------------------------------------

try:
    from pydantic import BaseModel as _BaseModel

    class SmileRequest(_BaseModel):
        """Request body containing a SMILES string."""

        smiles: str

except ImportError:
    # FastAPI / pydantic not installed; placeholder to avoid NameError at import
    SmileRequest = None  # type: ignore[assignment,misc]


def create_app() -> Any:
    """Create and configure the FastAPI application.

    Returns:
        Configured FastAPI application instance.

    Raises:
        ImportError: If FastAPI is not installed.
    """
    try:
        from fastapi import FastAPI, HTTPException, Query
    except ImportError as exc:
        raise ImportError(
            "FastAPI is required for the REST API. "
            "Install it with: pip install chemfuse[api]"
        ) from exc

    app = FastAPI(
        title="ChemFuse API",
        description="Multi-database cheminformatics REST API",
        version="0.2.0",
    )

    # ---------------------------------------------------------------------------
    # Helper: convert Compound to plain dict
    # ---------------------------------------------------------------------------

    def _compound_to_dict(compound: Any) -> dict[str, Any]:
        """Serialize a Compound model to a JSON-safe dict."""
        try:
            return compound.model_dump(exclude_none=True)
        except AttributeError:
            return dict(compound)

    # ---------------------------------------------------------------------------
    # Endpoints
    # ---------------------------------------------------------------------------

    @app.get("/health", summary="Health check")
    def health() -> dict[str, str]:
        """Return service health status and package version."""
        return {"status": "ok", "version": __version__}

    @app.get("/search", summary="Search compounds")
    def search(
        q: str = Query(..., description="Search query (name, SMILES, CID, formula, InChI)"),
        sources: str = Query(
            "pubchem",
            description="Comma-separated list of sources, e.g. 'pubchem,chembl'",
        ),
        query_type: str = Query(
            "name",
            description="Query type: name | smiles | cid | formula | inchi | inchikey",
        ),
        limit: int = Query(10, ge=1, le=200, description="Maximum results per source"),
    ) -> dict[str, Any]:
        """Search for compounds across one or more chemical databases."""
        import chemfuse

        source_list = [s.strip() for s in sources.split(",") if s.strip()]
        try:
            collection = chemfuse.search(
                q, sources=source_list, query_type=query_type, limit=limit
            )
        except Exception as exc:
            logger.exception("search endpoint error")
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        compounds = [_compound_to_dict(c) for c in collection.compounds]
        return {
            "compounds": compounds,
            "count": len(compounds),
            "sources": collection.sources,
        }

    @app.get("/compound/{identifier}", summary="Get compound by identifier")
    def get_compound(
        identifier: str,
        source: str = Query("pubchem", description="Source database"),
    ) -> dict[str, Any]:
        """Retrieve a specific compound by its database identifier."""
        import chemfuse

        try:
            compound = chemfuse.get(identifier, source=source)
        except Exception as exc:
            logger.exception("get_compound endpoint error for %s", identifier)
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        if compound is None:
            raise HTTPException(
                status_code=404,
                detail=f"Compound '{identifier}' not found in source '{source}'.",
            )
        return _compound_to_dict(compound)

    @app.get("/similar", summary="Find similar compounds")
    def find_similar(
        smiles: str = Query(..., description="Query SMILES string"),
        threshold: int = Query(90, ge=0, le=100, description="Tanimoto similarity threshold"),
        max_results: int = Query(10, ge=1, le=200, description="Maximum number of results"),
    ) -> dict[str, Any]:
        """Find compounds structurally similar to a given SMILES string."""
        import chemfuse

        if not smiles.strip():
            raise HTTPException(status_code=400, detail="SMILES string must not be empty.")
        try:
            collection = chemfuse.find_similar(
                smiles, threshold=threshold, max_results=max_results
            )
        except Exception as exc:
            logger.exception("find_similar endpoint error")
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        compounds = [_compound_to_dict(c) for c in collection.compounds]
        return {"compounds": compounds, "count": len(compounds)}

    @app.post("/admet", summary="Predict ADMET properties")
    def predict_admet(body: SmileRequest) -> dict[str, Any]:
        """Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties."""
        if not body.smiles.strip():
            raise HTTPException(status_code=400, detail="SMILES string must not be empty.")
        try:
            from chemfuse.compute.admet import predict_admet as _predict

            profile = _predict(body.smiles)
        except Exception as exc:
            logger.exception("admet endpoint error")
            detail = str(exc)
            if "Invalid SMILES" in detail or "invalid" in detail.lower():
                raise HTTPException(status_code=400, detail=detail) from exc
            raise HTTPException(status_code=500, detail=detail) from exc

        predictions: dict[str, Any] = {}
        for key, pred in profile.predictions.items():
            predictions[key] = {
                "value": pred.value,
                "unit": pred.unit,
                "confidence": pred.confidence,
                "method": pred.method,
                "category": pred.category,
            }

        return {
            "smiles": body.smiles,
            "predictions": predictions,
            "overall_score": profile.overall_score,
        }

    @app.post("/druglikeness", summary="Check drug-likeness")
    def check_druglikeness(body: SmileRequest) -> dict[str, Any]:
        """Evaluate drug-likeness using standard filters (Lipinski, Veber, Ghose, Egan, Muegge)."""
        if not body.smiles.strip():
            raise HTTPException(status_code=400, detail="SMILES string must not be empty.")
        try:
            from chemfuse.compute.druglikeness import check_drug_likeness as _check

            # Pass smiles=... so RDKit fills in properties automatically when available.
            # An empty properties dict lets the function use SMILES-derived values.
            dl = _check({}, smiles=body.smiles)
        except Exception as exc:
            logger.exception("druglikeness endpoint error")
            detail = str(exc)
            raise HTTPException(status_code=400, detail=detail) from exc

        def _filter_to_dict(f: Any) -> dict[str, Any]:
            return {
                "passed": f.pass_filter,
                "violations": f.violations,
                "details": f.details,
            }

        filters: dict[str, Any] = {}
        for name, f in {
            "lipinski": dl.lipinski,
            "veber": dl.veber,
            "ghose": dl.ghose,
            "egan": dl.egan,
            "muegge": dl.muegge,
        }.items():
            if f is not None:
                filters[name] = _filter_to_dict(f)

        if dl.pains is not None:
            filters["pains"] = {
                "passed": dl.pains.pass_filter,
                "violations": dl.pains.violations,
                "details": dl.pains.details,
            }
        if dl.qed is not None:
            filters["qed"] = dl.qed  # already a dict with "qed" and "classification" keys

        return {"filters": filters}

    @app.post("/descriptors", summary="Compute molecular descriptors")
    def compute_descriptors(body: SmileRequest) -> dict[str, Any]:
        """Compute 200+ RDKit molecular descriptors for a SMILES string."""
        if not body.smiles.strip():
            raise HTTPException(status_code=400, detail="SMILES string must not be empty.")
        try:
            from chemfuse.compute.descriptors import compute_descriptors as _compute

            descs = _compute(body.smiles)
        except ImportError as exc:
            raise HTTPException(
                status_code=400,
                detail="RDKit is required for descriptor computation. "
                "Install it with: pip install chemfuse[rdkit]",
            ) from exc
        except Exception as exc:
            logger.exception("descriptors endpoint error")
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        if not descs:
            raise HTTPException(
                status_code=400,
                detail=f"Could not parse SMILES string: {body.smiles!r}",
            )
        return {"smiles": body.smiles, "descriptors": descs}

    @app.get("/crossref", summary="Cross-reference identifiers")
    def crossref(
        identifier: str = Query(..., description="Compound identifier to look up"),
        type: str = Query(
            ...,
            description="Identifier type: cid | chembl_id | inchikey | smiles",
        ),
    ) -> dict[str, Any]:
        """Map a compound identifier across databases using UniChem."""
        valid_types = {"cid", "chembl_id", "inchikey", "smiles"}
        if type not in valid_types:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid identifier type '{type}'. Must be one of: {sorted(valid_types)}",
            )

        kwargs: dict[str, Any] = {}
        if type == "cid":
            try:
                kwargs["cid"] = int(identifier)
            except ValueError as exc:
                raise HTTPException(
                    status_code=400,
                    detail=f"CID must be an integer, got: {identifier!r}",
                ) from exc
        elif type == "chembl_id":
            kwargs["chembl_id"] = identifier
        elif type == "inchikey":
            kwargs["inchikey"] = identifier
        elif type == "smiles":
            kwargs["smiles"] = identifier

        try:
            import chemfuse

            mappings = chemfuse.map_identifiers(**kwargs)
        except Exception as exc:
            logger.exception("crossref endpoint error")
            raise HTTPException(status_code=500, detail=str(exc)) from exc

        return {"mappings": mappings}

    return app
