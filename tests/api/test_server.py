"""Tests for ChemFuse REST API server."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

try:
    from fastapi.testclient import TestClient

    from chemfuse.api.server import create_app

    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
INVALID_SMILES = "NOT_A_SMILES!!!"


def _make_compound(name: str = "aspirin", smiles: str = ASPIRIN_SMILES) -> MagicMock:
    """Return a minimal mock Compound that serialises via model_dump()."""
    c = MagicMock()
    c.model_dump.return_value = {
        "smiles": smiles,
        "name": name,
        "cid": 2244,
        "sources": ["pubchem"],
    }
    return c


def _make_collection(compounds: list | None = None, sources: list | None = None) -> MagicMock:
    col = MagicMock()
    col.compounds = compounds or [_make_compound()]
    col.sources = sources or ["pubchem"]
    return col


@pytest.fixture()
def client():
    """Return a TestClient for the ChemFuse FastAPI app."""
    return TestClient(create_app())


# ---------------------------------------------------------------------------
# Health endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestHealthEndpoint:
    def test_health_returns_ok(self, client: TestClient) -> None:
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json()["status"] == "ok"

    def test_health_contains_version(self, client: TestClient) -> None:
        response = client.get("/health")
        assert "version" in response.json()

    def test_health_version_is_string(self, client: TestClient) -> None:
        response = client.get("/health")
        assert isinstance(response.json()["version"], str)


# ---------------------------------------------------------------------------
# Search endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestSearchEndpoint:
    def test_search_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.search") as mock_search:
            mock_search.return_value = _make_collection()
            response = client.get("/search?q=aspirin")
        assert response.status_code == 200

    def test_search_response_structure(self, client: TestClient) -> None:
        with patch("chemfuse.search") as mock_search:
            mock_search.return_value = _make_collection()
            response = client.get("/search?q=aspirin")
        body = response.json()
        assert "compounds" in body
        assert "count" in body
        assert "sources" in body

    def test_search_count_matches_compounds(self, client: TestClient) -> None:
        with patch("chemfuse.search") as mock_search:
            mock_search.return_value = _make_collection()
            response = client.get("/search?q=aspirin")
        body = response.json()
        assert body["count"] == len(body["compounds"])

    def test_search_requires_q_param(self, client: TestClient) -> None:
        response = client.get("/search")
        assert response.status_code == 422  # Unprocessable Entity

    def test_search_passes_sources_param(self, client: TestClient) -> None:
        with patch("chemfuse.search") as mock_search:
            mock_search.return_value = _make_collection(sources=["pubchem", "chembl"])
            response = client.get("/search?q=aspirin&sources=pubchem,chembl")
        assert response.status_code == 200
        mock_search.assert_called_once()
        _, kwargs = mock_search.call_args
        assert kwargs["sources"] == ["pubchem", "chembl"]

    def test_search_passes_query_type(self, client: TestClient) -> None:
        with patch("chemfuse.search") as mock_search:
            mock_search.return_value = _make_collection()
            response = client.get("/search?q=CC(=O)O&query_type=smiles")
        assert response.status_code == 200
        _, kwargs = mock_search.call_args
        assert kwargs["query_type"] == "smiles"

    def test_search_handles_exception_as_500(self, client: TestClient) -> None:
        with patch("chemfuse.search", side_effect=RuntimeError("upstream failure")):
            response = client.get("/search?q=aspirin")
        assert response.status_code == 500


# ---------------------------------------------------------------------------
# Compound endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestCompoundEndpoint:
    def test_get_compound_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.get") as mock_get:
            mock_get.return_value = _make_compound()
            response = client.get("/compound/2244")
        assert response.status_code == 200

    def test_get_compound_not_found_returns_404(self, client: TestClient) -> None:
        with patch("chemfuse.get", return_value=None):
            response = client.get("/compound/99999999")
        assert response.status_code == 404

    def test_get_compound_contains_smiles(self, client: TestClient) -> None:
        with patch("chemfuse.get") as mock_get:
            mock_get.return_value = _make_compound()
            response = client.get("/compound/2244")
        body = response.json()
        assert "smiles" in body

    def test_get_compound_passes_source(self, client: TestClient) -> None:
        with patch("chemfuse.get") as mock_get:
            mock_get.return_value = _make_compound()
            response = client.get("/compound/CHEMBL25?source=chembl")
        assert response.status_code == 200
        _, kwargs = mock_get.call_args
        assert kwargs["source"] == "chembl"


# ---------------------------------------------------------------------------
# Similar endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestSimilarEndpoint:
    def test_similar_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.find_similar") as mock_sim:
            mock_sim.return_value = _make_collection()
            response = client.get(f"/similar?smiles={ASPIRIN_SMILES}")
        assert response.status_code == 200

    def test_similar_response_structure(self, client: TestClient) -> None:
        with patch("chemfuse.find_similar") as mock_sim:
            mock_sim.return_value = _make_collection()
            response = client.get(f"/similar?smiles={ASPIRIN_SMILES}")
        body = response.json()
        assert "compounds" in body
        assert "count" in body

    def test_similar_requires_smiles(self, client: TestClient) -> None:
        response = client.get("/similar")
        assert response.status_code == 422

    def test_similar_empty_smiles_returns_400(self, client: TestClient) -> None:
        response = client.get("/similar?smiles=%20")  # URL-encoded space
        assert response.status_code == 400

    def test_similar_passes_threshold(self, client: TestClient) -> None:
        with patch("chemfuse.find_similar") as mock_sim:
            mock_sim.return_value = _make_collection()
            response = client.get(f"/similar?smiles={ASPIRIN_SMILES}&threshold=75")
        assert response.status_code == 200
        _, kwargs = mock_sim.call_args
        assert kwargs["threshold"] == 75


# ---------------------------------------------------------------------------
# ADMET endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestAdmetEndpoint:
    def _make_profile(self) -> MagicMock:
        pred = MagicMock()
        pred.value = 0.8
        pred.unit = None
        pred.confidence = 0.7
        pred.method = "rule-based"
        pred.category = "high"

        profile = MagicMock()
        profile.predictions = {"solubility": pred}
        profile.overall_score = 0.75
        return profile

    def test_admet_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.compute.admet.predict_admet") as mock_pred:
            mock_pred.return_value = self._make_profile()
            response = client.post("/admet", json={"smiles": ASPIRIN_SMILES})
        assert response.status_code == 200

    def test_admet_response_structure(self, client: TestClient) -> None:
        with patch("chemfuse.compute.admet.predict_admet") as mock_pred:
            mock_pred.return_value = self._make_profile()
            response = client.post("/admet", json={"smiles": ASPIRIN_SMILES})
        body = response.json()
        assert "smiles" in body
        assert "predictions" in body
        assert "overall_score" in body

    def test_admet_empty_smiles_returns_400(self, client: TestClient) -> None:
        response = client.post("/admet", json={"smiles": ""})
        assert response.status_code == 400

    def test_admet_missing_body_returns_422(self, client: TestClient) -> None:
        response = client.post("/admet", json={})
        assert response.status_code == 422

    def test_admet_real_aspirin(self, client: TestClient) -> None:
        """Integration-style test using real compute module (RDKit optional)."""
        response = client.post("/admet", json={"smiles": ASPIRIN_SMILES})
        # Must succeed or return 400 (invalid SMILES) / 500 (missing RDKit)
        assert response.status_code in (200, 400, 500)
        if response.status_code == 200:
            body = response.json()
            assert body["smiles"] == ASPIRIN_SMILES
            assert isinstance(body["predictions"], dict)


# ---------------------------------------------------------------------------
# Drug-likeness endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestDruglikenessEndpoint:
    def _make_druglikeness(self) -> MagicMock:
        def _filter(passed: bool = True) -> MagicMock:
            f = MagicMock()
            f.pass_filter = passed
            f.violations = []
            f.details = {}
            return f

        dl = MagicMock()
        dl.lipinski = _filter()
        dl.veber = _filter()
        dl.ghose = _filter()
        dl.egan = _filter()
        dl.muegge = _filter()
        dl.pains = None
        dl.qed = None
        return dl

    def test_druglikeness_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.compute.druglikeness.check_drug_likeness") as mock_dl:
            mock_dl.return_value = self._make_druglikeness()
            response = client.post("/druglikeness", json={"smiles": ASPIRIN_SMILES})
        assert response.status_code == 200

    def test_druglikeness_response_structure(self, client: TestClient) -> None:
        with patch("chemfuse.compute.druglikeness.check_drug_likeness") as mock_dl:
            mock_dl.return_value = self._make_druglikeness()
            response = client.post("/druglikeness", json={"smiles": ASPIRIN_SMILES})
        body = response.json()
        assert "filters" in body

    def test_druglikeness_filters_contain_lipinski(self, client: TestClient) -> None:
        with patch("chemfuse.compute.druglikeness.check_drug_likeness") as mock_dl:
            mock_dl.return_value = self._make_druglikeness()
            response = client.post("/druglikeness", json={"smiles": ASPIRIN_SMILES})
        filters = response.json()["filters"]
        assert "lipinski" in filters

    def test_druglikeness_empty_smiles_returns_400(self, client: TestClient) -> None:
        response = client.post("/druglikeness", json={"smiles": ""})
        assert response.status_code == 400

    def test_druglikeness_real_aspirin(self, client: TestClient) -> None:
        """Integration-style test using real compute module."""
        response = client.post("/druglikeness", json={"smiles": ASPIRIN_SMILES})
        assert response.status_code in (200, 400, 500)
        if response.status_code == 200:
            body = response.json()
            assert "filters" in body
            assert "lipinski" in body["filters"]


# ---------------------------------------------------------------------------
# Descriptors endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestDescriptorsEndpoint:
    def test_descriptors_empty_smiles_returns_400(self, client: TestClient) -> None:
        response = client.post("/descriptors", json={"smiles": ""})
        assert response.status_code == 400

    def test_descriptors_missing_smiles_returns_422(self, client: TestClient) -> None:
        response = client.post("/descriptors", json={})
        assert response.status_code == 422

    def test_descriptors_mocked_success(self, client: TestClient) -> None:
        mock_descs = {"MolWt": 180.16, "MolLogP": 1.19, "NumHDonors": 1}
        with patch("chemfuse.compute.descriptors.compute_descriptors", return_value=mock_descs):
            response = client.post("/descriptors", json={"smiles": ASPIRIN_SMILES})
        assert response.status_code == 200
        body = response.json()
        assert body["smiles"] == ASPIRIN_SMILES
        assert "descriptors" in body
        assert body["descriptors"]["MolWt"] == pytest.approx(180.16)

    def test_descriptors_empty_result_returns_400(self, client: TestClient) -> None:
        """An invalid SMILES that returns {} from compute_descriptors gives 400."""
        with patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}):
            response = client.post("/descriptors", json={"smiles": INVALID_SMILES})
        assert response.status_code == 400

    def test_descriptors_rdkit_not_installed_returns_400(self, client: TestClient) -> None:
        with patch(
            "chemfuse.compute.descriptors.compute_descriptors",
            side_effect=ImportError("rdkit not found"),
        ):
            response = client.post("/descriptors", json={"smiles": ASPIRIN_SMILES})
        assert response.status_code == 400


# ---------------------------------------------------------------------------
# Cross-reference endpoint
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not FASTAPI_AVAILABLE, reason="FastAPI not installed")
class TestCrossrefEndpoint:
    def test_crossref_by_inchikey_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.map_identifiers", return_value={"pubchem": "2244"}):
            response = client.get(
                "/crossref?identifier=BSYNRYMUTXBXSQ-UHFFFAOYSA-N&type=inchikey"
            )
        assert response.status_code == 200
        body = response.json()
        assert "mappings" in body

    def test_crossref_by_cid_returns_200(self, client: TestClient) -> None:
        with patch("chemfuse.map_identifiers", return_value={"pubchem": "2244"}):
            response = client.get("/crossref?identifier=2244&type=cid")
        assert response.status_code == 200

    def test_crossref_invalid_type_returns_400(self, client: TestClient) -> None:
        response = client.get("/crossref?identifier=foo&type=badtype")
        assert response.status_code == 400

    def test_crossref_cid_non_integer_returns_400(self, client: TestClient) -> None:
        response = client.get("/crossref?identifier=notanumber&type=cid")
        assert response.status_code == 400

    def test_crossref_missing_params_returns_422(self, client: TestClient) -> None:
        response = client.get("/crossref")
        assert response.status_code == 422

    def test_crossref_mappings_is_dict(self, client: TestClient) -> None:
        with patch(
            "chemfuse.map_identifiers",
            return_value={"pubchem": "2244", "chembl": "CHEMBL25"},
        ):
            response = client.get("/crossref?identifier=2244&type=cid")
        body = response.json()
        assert isinstance(body["mappings"], dict)
