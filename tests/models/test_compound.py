"""Tests for chemfuse.models.compound."""

from __future__ import annotations

import pytest

from chemfuse.models.compound import Compound, CompoundProperties


class TestCompoundProperties:
    def test_all_none_by_default(self):
        props = CompoundProperties()
        assert props.molecular_weight is None
        assert props.xlogp is None
        assert props.tpsa is None

    def test_construct_with_values(self):
        props = CompoundProperties(
            molecular_weight=180.16,
            xlogp=1.2,
            tpsa=63.6,
            hbd_count=1,
            hba_count=4,
        )
        assert props.molecular_weight == 180.16
        assert props.xlogp == 1.2
        assert props.tpsa == 63.6
        assert props.hbd_count == 1
        assert props.hba_count == 4

    def test_invalid_type_raises(self):
        from pydantic import ValidationError as PydanticValidationError

        with pytest.raises(PydanticValidationError):
            CompoundProperties(molecular_weight="not-a-number")  # type: ignore[arg-type]


class TestCompoundCreation:
    def test_minimal_compound(self):
        c = Compound(smiles="C")
        assert c.smiles == "C"
        assert c.cid is None
        assert c.name is None
        assert c.sources == []

    def test_full_compound(self):
        props = CompoundProperties(molecular_weight=180.16)
        c = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchi="InChI=1S/C9H8O4",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            name="aspirin",
            formula="C9H8O4",
            synonyms=["aspirin", "acetylsalicylic acid"],
            sources=["pubchem"],
            properties=props,
        )
        assert c.cid == 2244
        assert c.name == "aspirin"
        assert c.formula == "C9H8O4"
        assert "pubchem" in c.sources
        assert c.properties.molecular_weight == 180.16

    def test_default_properties_empty(self):
        c = Compound(smiles="C")
        assert isinstance(c.properties, CompoundProperties)
        assert c.properties.molecular_weight is None

    def test_repr_with_name(self, aspirin: Compound):
        r = repr(aspirin)
        assert "aspirin" in r

    def test_repr_with_smiles_only(self):
        c = Compound(smiles="CC(=O)O")
        r = repr(c)
        assert "Compound" in r


class TestCompoundToDict:
    def test_to_dict_basic(self, aspirin: Compound):
        d = aspirin.to_dict()
        assert isinstance(d, dict)
        assert "smiles" in d
        assert "name" in d

    def test_to_dict_flattens_properties(self, aspirin: Compound):
        d = aspirin.to_dict()
        # Properties should be flattened into top-level keys
        assert "molecular_weight" in d
        assert "properties" not in d

    def test_to_dict_skips_none_properties(self):
        c = Compound(
            smiles="C",
            properties=CompoundProperties(molecular_weight=180.0),
        )
        d = c.to_dict()
        assert "molecular_weight" in d
        # None property values should not appear
        assert "xlogp" not in d or d.get("xlogp") is not None  # Only if set

    def test_to_dict_includes_sources(self, aspirin: Compound):
        d = aspirin.to_dict()
        assert "sources" in d
        assert "pubchem" in d["sources"]


class TestCompoundToSeries:
    def test_to_series_returns_series(self, aspirin: Compound):
        import pandas as pd
        s = aspirin.to_series()
        assert isinstance(s, pd.Series)

    def test_to_series_has_name(self, aspirin: Compound):
        s = aspirin.to_series()
        assert s["name"] == "aspirin"


class TestCompoundMerge:
    def test_merge_by_inchikey(self):
        c1 = Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            name="aspirin",
            sources=["pubchem"],
        )
        c2 = Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            formula="C9H8O4",
            sources=["chembl"],
            properties=CompoundProperties(molecular_weight=180.16),
        )
        merged = c1.merge(c2)
        assert merged.formula == "C9H8O4"
        assert merged.name == "aspirin"
        assert "pubchem" in merged.sources
        assert "chembl" in merged.sources
        assert merged.properties.molecular_weight == 180.16

    def test_merge_by_cid(self):
        c1 = Compound(cid=2244, smiles="CC(=O)Oc1ccccc1C(=O)O", name="aspirin", sources=["pubchem"])
        c2 = Compound(cid=2244, smiles="CC(=O)Oc1ccccc1C(=O)O", formula="C9H8O4", sources=["other"])
        merged = c1.merge(c2)
        assert merged.formula == "C9H8O4"
        assert merged.name == "aspirin"

    def test_merge_by_smiles(self):
        c1 = Compound(smiles="CCO", name="ethanol", sources=["src1"])
        c2 = Compound(smiles="CCO", formula="C2H6O", sources=["src2"])
        merged = c1.merge(c2)
        assert merged.formula == "C2H6O"
        assert merged.name == "ethanol"

    def test_merge_does_not_overwrite_existing_fields(self):
        c1 = Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            name="aspirin",
            sources=["pubchem"],
        )
        c2 = Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            name="acetylsalicylic acid",
            sources=["chembl"],
        )
        merged = c1.merge(c2)
        # c1's name should take priority
        assert merged.name == "aspirin"

    def test_merge_synonyms_combined(self):
        c1 = Compound(
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            synonyms=["ethanol"],
            sources=["src1"],
        )
        c2 = Compound(
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            synonyms=["alcohol", "ethyl alcohol"],
            sources=["src2"],
        )
        merged = c1.merge(c2)
        assert "ethanol" in merged.synonyms
        assert "alcohol" in merged.synonyms

    def test_merge_different_molecules_raises(self):
        c1 = Compound(
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            sources=["pubchem"],
        )
        c2 = Compound(
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            sources=["pubchem"],
        )
        with pytest.raises(ValueError):
            c1.merge(c2)

    def test_merge_fills_missing_properties(self):
        c1 = Compound(
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            sources=["src1"],
        )
        c2 = Compound(
            smiles="CCO",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            sources=["src2"],
            properties=CompoundProperties(molecular_weight=46.07, xlogp=-0.3),
        )
        merged = c1.merge(c2)
        assert merged.properties.molecular_weight == 46.07
        assert merged.properties.xlogp == -0.3


class TestCompoundCFE08Fields:
    """Tests for CF-E08: clinical/regulatory and cross-reference fields on Compound."""

    def test_clinical_fields_default_none(self) -> None:
        """max_phase, molecule_type, is_withdrawn, black_box_warning default to None."""
        c = Compound(smiles="CCO")
        assert c.max_phase is None
        assert c.molecule_type is None
        assert c.is_withdrawn is None
        assert c.black_box_warning is None

    def test_cross_reference_fields_default_none(self) -> None:
        """drugbank_id, kegg_id, chebi_id default to None."""
        c = Compound(smiles="CCO")
        assert c.drugbank_id is None
        assert c.kegg_id is None
        assert c.chebi_id is None

    def test_max_phase_stored(self) -> None:
        """max_phase=4 (approved) is stored correctly."""
        c = Compound(smiles="CC(=O)Oc1ccccc1C(O)=O", max_phase=4)
        assert c.max_phase == 4

    def test_max_phase_zero_is_not_none(self) -> None:
        """max_phase=0 (preclinical) is stored as 0, not treated as falsy/None."""
        c = Compound(smiles="CCO", max_phase=0)
        assert c.max_phase == 0
        assert c.max_phase is not None

    def test_molecule_type_stored(self) -> None:
        """molecule_type is stored correctly."""
        c = Compound(smiles="CCO", molecule_type="Small molecule")
        assert c.molecule_type == "Small molecule"

    def test_is_withdrawn_flag(self) -> None:
        """is_withdrawn bool flag is stored."""
        c = Compound(smiles="CCO", is_withdrawn=True)
        assert c.is_withdrawn is True

    def test_black_box_warning_flag(self) -> None:
        """black_box_warning bool flag is stored."""
        c = Compound(smiles="CCO", black_box_warning=True)
        assert c.black_box_warning is True

    def test_drugbank_id_stored(self) -> None:
        """drugbank_id is stored correctly."""
        c = Compound(smiles="CC(=O)Oc1ccccc1C(O)=O", drugbank_id="DB00945")
        assert c.drugbank_id == "DB00945"

    def test_kegg_id_stored(self) -> None:
        """kegg_id is stored correctly."""
        c = Compound(smiles="CC(=O)Oc1ccccc1C(O)=O", kegg_id="D00109")
        assert c.kegg_id == "D00109"

    def test_chebi_id_stored(self) -> None:
        """chebi_id is stored correctly."""
        c = Compound(smiles="CC(=O)Oc1ccccc1C(O)=O", chebi_id="CHEBI:15365")
        assert c.chebi_id == "CHEBI:15365"

    def test_cf_e08_fields_in_model_dump(self) -> None:
        """All CF-E08 fields appear in model_dump output."""
        c = Compound(
            smiles="CCO",
            max_phase=4,
            molecule_type="Small molecule",
            is_withdrawn=False,
            black_box_warning=False,
            drugbank_id="DB00945",
            kegg_id="D00109",
            chebi_id="CHEBI:15365",
        )
        d = c.model_dump()
        assert d["max_phase"] == 4
        assert d["molecule_type"] == "Small molecule"
        assert d["is_withdrawn"] is False
        assert d["black_box_warning"] is False
        assert d["drugbank_id"] == "DB00945"
        assert d["kegg_id"] == "D00109"
        assert d["chebi_id"] == "CHEBI:15365"
