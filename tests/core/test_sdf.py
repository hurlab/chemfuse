"""Tests for chemfuse.core.sdf (CF-E09: SDF Import + Enriched SDF Export)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from chemfuse.core.sdf import read_sdf, write_sdf
from chemfuse.models.compound import Compound, CompoundProperties

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_sdf(tmp_path: Path) -> str:
    """Create a small SDF file with two molecules for testing."""
    from rdkit import Chem

    writer = Chem.SDWriter(str(tmp_path / "test.sdf"))

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol.SetProp("_Name", "aspirin")
    mol.SetProp("CID", "2244")
    writer.write(mol)

    mol2 = Chem.MolFromSmiles("CCO")
    mol2.SetProp("_Name", "ethanol")
    writer.write(mol2)

    writer.close()
    return str(tmp_path / "test.sdf")


@pytest.fixture
def three_compound_sdf(tmp_path: Path) -> str:
    """SDF with three compounds including extra SD properties."""
    from rdkit import Chem

    writer = Chem.SDWriter(str(tmp_path / "three.sdf"))

    smiles_list = [
        ("caffeine", "Cn1c(=O)c2c(ncn2C)n(c1=O)C"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ]
    for name, smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        mol.SetProp("_Name", name)
        writer.write(mol)

    writer.close()
    return str(tmp_path / "three.sdf")


@pytest.fixture
def enriched_compounds() -> list[Compound]:
    """Three Compound objects with identifiers and properties set."""
    return [
        Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            name="aspirin",
            cid=2244,
            chembl_id="CHEMBL25",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            sources=["pubchem", "chembl"],
            properties=CompoundProperties(
                molecular_weight=180.16,
                xlogp=1.2,
                tpsa=63.6,
                hbd_count=1,
                hba_count=4,
            ),
        ),
        Compound(
            smiles="CCO",
            name="ethanol",
            cid=702,
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=46.07),
        ),
    ]


# ---------------------------------------------------------------------------
# read_sdf tests
# ---------------------------------------------------------------------------

class TestReadSdf:
    def test_reads_three_compounds(self, three_compound_sdf: str) -> None:
        compounds = read_sdf(three_compound_sdf)
        assert len(compounds) == 3

    def test_smiles_extracted(self, sample_sdf: str) -> None:
        compounds = read_sdf(sample_sdf)
        smiles_set = {c.smiles for c in compounds}
        assert len(smiles_set) == 2
        assert all(s for s in smiles_set)

    def test_name_extracted(self, sample_sdf: str) -> None:
        compounds = read_sdf(sample_sdf)
        names = {c.name for c in compounds}
        assert "aspirin" in names
        assert "ethanol" in names

    def test_cid_sd_tag_mapped(self, sample_sdf: str) -> None:
        compounds = read_sdf(sample_sdf)
        aspirin = next(c for c in compounds if c.name == "aspirin")
        assert aspirin.cid == 2244

    def test_sources_set_to_sdf(self, sample_sdf: str) -> None:
        compounds = read_sdf(sample_sdf)
        for compound in compounds:
            assert "sdf" in compound.sources

    def test_empty_sdf_returns_empty_list(self, tmp_path: Path) -> None:
        # RDKit requires at least a valid (but record-less) SDF: just a blank line
        empty_path = tmp_path / "empty.sdf"
        empty_path.write_text("\n")

        compounds = read_sdf(str(empty_path))
        assert compounds == []

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="SDF file not found"):
            read_sdf(str(tmp_path / "nonexistent.sdf"))

    def test_optional_dependency_error_when_rdkit_missing(self, sample_sdf: str) -> None:
        import builtins

        real_import = builtins.__import__

        def mock_import(name: str, *args, **kwargs):
            if name == "rdkit" or name.startswith("rdkit."):
                raise ImportError("Mocked rdkit absence")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=mock_import):
            from chemfuse.core.exceptions import OptionalDependencyError

            with pytest.raises(OptionalDependencyError):
                read_sdf(sample_sdf)

    def test_chembl_id_sd_tag_mapped(self, tmp_path: Path) -> None:
        from rdkit import Chem

        path = str(tmp_path / "with_chembl.sdf")
        writer = Chem.SDWriter(path)
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol.SetProp("_Name", "benzene")
        mol.SetProp("ChEMBL_ID", "CHEMBL277500")
        writer.write(mol)
        writer.close()

        compounds = read_sdf(path)
        assert compounds[0].chembl_id == "CHEMBL277500"

    def test_mw_sd_tag_mapped_to_properties(self, tmp_path: Path) -> None:
        from rdkit import Chem

        path = str(tmp_path / "with_mw.sdf")
        writer = Chem.SDWriter(path)
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        mol.SetProp("_Name", "aspirin")
        mol.SetProp("MW", "180.16")
        writer.write(mol)
        writer.close()

        compounds = read_sdf(path)
        assert compounds[0].properties.molecular_weight == pytest.approx(180.16)


# ---------------------------------------------------------------------------
# write_sdf tests
# ---------------------------------------------------------------------------

class TestWriteSdf:
    def test_creates_file(self, enriched_compounds: list[Compound], tmp_path: Path) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)
        assert Path(out).exists()

    def test_returns_path_string(self, enriched_compounds: list[Compound], tmp_path: Path) -> None:
        out = str(tmp_path / "out.sdf")
        result = write_sdf(enriched_compounds, out)
        assert isinstance(result, str)
        assert Path(result).exists()

    def test_writes_correct_molecule_count(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mols = [m for m in supplier if m is not None]
        assert len(mols) == len(enriched_compounds)

    def test_sd_tags_include_cid(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol is not None
        assert mol.HasProp("CID")
        assert mol.GetProp("CID") == "2244"

    def test_sd_tags_include_chembl_id(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol.HasProp("ChEMBL_ID")
        assert mol.GetProp("ChEMBL_ID") == "CHEMBL25"

    def test_sd_tags_include_molecular_weight(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol.HasProp("MW")
        assert float(mol.GetProp("MW")) == pytest.approx(180.16)

    def test_sd_tags_include_sources(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "out.sdf")
        write_sdf(enriched_compounds, out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol.HasProp("sources")
        sources = mol.GetProp("sources").split(",")
        assert "pubchem" in sources
        assert "chembl" in sources

    def test_include_properties_false_produces_minimal_sdf(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "minimal.sdf")
        write_sdf(enriched_compounds, out, include_properties=False)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol is not None
        # No data properties should be present (only _Name is set by SDWriter)
        prop_names = [p for p in mol.GetPropNames() if not p.startswith("_")]
        assert prop_names == []

    def test_descriptors_written_with_desc_prefix(self, tmp_path: Path) -> None:
        compound = Compound(
            smiles="c1ccccc1",
            name="benzene",
            sources=["test"],
            descriptors={"MolLogP": 2.1, "TPSA": 0.0},
        )
        out = str(tmp_path / "descriptors.sdf")
        write_sdf([compound], out)

        from rdkit import Chem

        supplier = Chem.SDMolSupplier(out)
        mol = next(iter(supplier))
        assert mol.HasProp("desc_MolLogP")
        assert float(mol.GetProp("desc_MolLogP")) == pytest.approx(2.1)

    def test_optional_dependency_error_when_rdkit_missing(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        import builtins

        real_import = builtins.__import__

        def mock_import(name: str, *args, **kwargs):
            if name == "rdkit" or name.startswith("rdkit."):
                raise ImportError("Mocked rdkit absence")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=mock_import):
            from chemfuse.core.exceptions import OptionalDependencyError

            with pytest.raises(OptionalDependencyError):
                write_sdf(enriched_compounds, str(tmp_path / "x.sdf"))


# ---------------------------------------------------------------------------
# Round-trip test
# ---------------------------------------------------------------------------

class TestRoundTrip:
    def test_write_then_read_preserves_smiles(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "roundtrip.sdf")
        write_sdf(enriched_compounds, out)
        read_back = read_sdf(out)

        assert len(read_back) == len(enriched_compounds)
        written_smiles = {c.smiles for c in enriched_compounds}
        read_smiles = {c.smiles for c in read_back}
        assert written_smiles == read_smiles

    def test_write_then_read_preserves_name(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "roundtrip_name.sdf")
        write_sdf(enriched_compounds, out)
        read_back = read_sdf(out)

        written_names = {c.name for c in enriched_compounds if c.name}
        read_names = {c.name for c in read_back if c.name}
        assert written_names == read_names

    def test_write_then_read_preserves_cid(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "roundtrip_cid.sdf")
        write_sdf(enriched_compounds, out)
        read_back = read_sdf(out)

        written_cids = {c.cid for c in enriched_compounds if c.cid is not None}
        read_cids = {c.cid for c in read_back if c.cid is not None}
        assert written_cids == read_cids

    def test_write_then_read_preserves_molecular_weight(
        self, enriched_compounds: list[Compound], tmp_path: Path
    ) -> None:
        out = str(tmp_path / "roundtrip_mw.sdf")
        write_sdf(enriched_compounds, out)
        read_back = read_sdf(out)

        for original, recovered in zip(enriched_compounds, read_back, strict=True):
            if original.properties.molecular_weight is not None:
                assert recovered.properties.molecular_weight == pytest.approx(
                    original.properties.molecular_weight, rel=1e-4
                )
