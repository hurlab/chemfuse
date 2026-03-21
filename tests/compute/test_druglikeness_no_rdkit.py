"""Tests for drug-likeness filters operating without RDKit.

Simulates the scenario where only PubChem-fetched properties are available
(no RDKit installed). Verifies that all five filters produce the same
pass/fail decisions as when RDKit is available.
"""

from __future__ import annotations

ASPIRIN_PROPS = {
    "molecular_weight": 180.16,
    "xlogp": 1.2,
    "tpsa": 63.6,
    "hbd_count": 1,
    "hba_count": 4,
    "rotatable_bonds": 3,
    "heavy_atom_count": 13,
}

CAFFEINE_PROPS = {
    "molecular_weight": 194.19,
    "xlogp": -0.07,
    "tpsa": 58.4,
    "hbd_count": 0,
    "hba_count": 6,
    "rotatable_bonds": 0,
    "heavy_atom_count": 14,
}

CYCLOSPORINE_PROPS = {
    "molecular_weight": 1202.61,
    "xlogp": 2.9,
    "tpsa": 278.8,
    "hbd_count": 5,
    "hba_count": 23,
    "rotatable_bonds": 15,
    "heavy_atom_count": 88,
}


def _disable_rdkit(module):
    """Context manager that disables RDKIT_AVAILABLE in a module."""
    import contextlib

    @contextlib.contextmanager
    def _ctx():
        original = module._RDKIT_AVAILABLE
        module._RDKIT_AVAILABLE = False
        try:
            yield
        finally:
            module._RDKIT_AVAILABLE = original

    return _ctx()


class TestLipinskiNoRDKit:
    def test_aspirin_passes(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import lipinski_filter
        with _disable_rdkit(dl_mod):
            result = lipinski_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_caffeine_passes(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import lipinski_filter
        with _disable_rdkit(dl_mod):
            result = lipinski_filter(CAFFEINE_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import lipinski_filter
        with _disable_rdkit(dl_mod):
            result = lipinski_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False


class TestVeberNoRDKit:
    def test_aspirin_passes(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import veber_filter
        with _disable_rdkit(dl_mod):
            result = veber_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import veber_filter
        with _disable_rdkit(dl_mod):
            result = veber_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False


class TestGhoseNoRDKit:
    def test_aspirin_passes(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import ghose_filter
        with _disable_rdkit(dl_mod):
            result = ghose_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import ghose_filter
        with _disable_rdkit(dl_mod):
            result = ghose_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False


class TestEganNoRDKit:
    def test_aspirin_passes(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import egan_filter
        with _disable_rdkit(dl_mod):
            result = egan_filter(ASPIRIN_PROPS)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import egan_filter
        with _disable_rdkit(dl_mod):
            result = egan_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False


class TestMueggeNoRDKit:
    def test_aspirin_no_rdkit(self):
        """Without RDKit, Muegge needs ring and heteroatom data from properties."""
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import muegge_filter
        # Aspirin minimal props that Muegge can evaluate
        props = {
            "molecular_weight": 180.16,
            "xlogp": 1.2,
            "tpsa": 63.6,
            "hbd_count": 1,
            "hba_count": 4,
            "rotatable_bonds": 3,
            "ring_count": 1,
            "heteroatom_count": 4,
            "carbon_count": 9,
        }
        with _disable_rdkit(dl_mod):
            result = muegge_filter(props)
        assert result.pass_filter is True

    def test_cyclosporine_fails(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import muegge_filter
        with _disable_rdkit(dl_mod):
            result = muegge_filter(CYCLOSPORINE_PROPS)
        assert result.pass_filter is False


class TestCheckDrugLikenessNoRDKit:
    def test_all_filters_work_without_rdkit(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import check_drug_likeness
        from chemfuse.models.prediction import DrugLikeness
        with _disable_rdkit(dl_mod):
            result = check_drug_likeness(ASPIRIN_PROPS)
        assert isinstance(result, DrugLikeness)

    def test_no_smiles_no_rdkit_uses_provided_props(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import check_drug_likeness
        with _disable_rdkit(dl_mod):
            result = check_drug_likeness(ASPIRIN_PROPS, smiles=None)
        # Should produce valid filter results using only provided data
        assert result.lipinski.pass_filter is True
        assert result.veber.pass_filter is True

    def test_rdkit_and_no_rdkit_agree_on_lipinski_aspirin(self):
        """RDKit and no-RDKit paths should agree on Lipinski for aspirin."""
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import lipinski_filter

        result_with = lipinski_filter(ASPIRIN_PROPS, smiles="CC(=O)Oc1ccccc1C(=O)O")
        with _disable_rdkit(dl_mod):
            result_without = lipinski_filter(ASPIRIN_PROPS)

        assert result_with.pass_filter == result_without.pass_filter

    def test_rdkit_and_no_rdkit_agree_on_lipinski_cyclosporine(self):
        import chemfuse.compute.druglikeness as dl_mod
        from chemfuse.compute.druglikeness import lipinski_filter

        result_with = lipinski_filter(CYCLOSPORINE_PROPS)
        with _disable_rdkit(dl_mod):
            result_without = lipinski_filter(CYCLOSPORINE_PROPS)

        assert result_with.pass_filter == result_without.pass_filter
