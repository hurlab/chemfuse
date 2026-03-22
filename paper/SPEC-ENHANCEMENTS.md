# ChemFuse Enhancement SPECs (CF-E01 through CF-E10)

## CF-E01: Murcko Scaffold Decomposition + Frequency Analysis

**Goal**: Add scaffold awareness to ChemFuse — extract Bemis-Murcko scaffolds, compute generic scaffolds, and provide scaffold frequency analysis on collections.

**Files to create/modify**:
- NEW: `src/chemfuse/analyze/scaffolds.py`
- MODIFY: `src/chemfuse/analyze/__init__.py` (export new functions)
- MODIFY: `src/chemfuse/models/compound.py` (add `scaffold` and `generic_scaffold` fields)
- MODIFY: `src/chemfuse/models/collection.py` (add `scaffold_frequency()` and `group_by_scaffold()` methods)
- NEW: `tests/analyze/test_scaffolds.py`

**API Design**:
```python
# analyze/scaffolds.py
def murcko_scaffold(smiles: str) -> str | None
def generic_scaffold(smiles: str) -> str | None
def scaffold_frequency(smiles_list: list[str]) -> dict[str, int]
def group_by_scaffold(compounds: list[Compound]) -> dict[str, list[Compound]]

# On CompoundCollection:
collection.scaffold_frequency() -> pd.DataFrame  # scaffold_smiles, count, percentage
collection.group_by_scaffold() -> dict[str, CompoundCollection]
```

**RDKit modules**: `rdkit.Chem.Scaffolds.MurckoScaffold` (GetScaffoldForMol, MakeScaffoldGeneric)

**Tests**: scaffold extraction for aspirin, caffeine, ibuprofen; generic scaffold strips heteroatoms; frequency counts correct; empty/invalid SMILES handled; collection methods work.

---

## CF-E02: Bioactivity Unit Normalization + pIC50 Conversion

**Goal**: Normalize all bioactivity values to nM and provide pIC50/pKi computed properties. Fix the silent data quality defect where mixed-unit values corrupt downstream analysis.

**Files to create/modify**:
- MODIFY: `src/chemfuse/models/bioactivity.py` (add `value_nm`, `pic50`, normalization logic)
- MODIFY: `src/chemfuse/sources/chembl.py` (normalize at parse time)
- MODIFY: `src/chemfuse/sources/bindingdb.py` (normalize at parse time)
- MODIFY: `src/chemfuse/models/compound.py` (add `best_activity(target, type)` method)
- MODIFY: `tests/models/test_bioactivity.py`
- MODIFY: `tests/sources/test_chembl.py`

**API Design**:
```python
# bioactivity.py additions
class Bioactivity:
    value_nm: float | None = None      # normalized to nM
    pic50: float | None = None          # -log10(value_nm * 1e-9), only for IC50
    pki: float | None = None            # -log10(value_nm * 1e-9), only for Ki
    is_normalized: bool = False

def normalize_to_nm(value: float, units: str) -> float | None
```

**Unit conversion table**: M->1e9, mM->1e6, uM/µM->1e3, nM->1, pM->1e-3. Handle `ug/mL` via MW if available.

**Tests**: 100 nM stays 100; 0.1 µM becomes 100; pIC50 of 100 nM = 7.0; None units return None; mixed-unit collection sorts correctly.

---

## CF-E03: Compound Standardization (Salt Stripping + Tautomer Canonicalization)

**Goal**: Add structure standardization to ensure correct deduplication and consistent descriptor computation.

**Files to create/modify**:
- NEW: `src/chemfuse/compute/standardization.py`
- MODIFY: `src/chemfuse/compute/__init__.py` (export)
- MODIFY: `src/chemfuse/models/compound.py` (add `to_mol()` method, `canonical_smiles` field)
- NEW: `tests/compute/test_standardization.py`

**API Design**:
```python
# standardization.py
def standardize_mol(smiles: str, strip_salts: bool = True, canonicalize_tautomers: bool = True) -> str | None
def strip_salts(smiles: str) -> str | None
def canonical_tautomer(smiles: str) -> str | None

# On Compound:
def to_mol(self) -> Any  # returns rdkit.Chem.Mol or None
```

**RDKit modules**: `rdkit.Chem.MolStandardize.rdMolStandardize` (Uncharger, LargestFragmentChooser, TautomerEnumerator, TautomerCanonicalizer)

**Tests**: aspirin sodium salt -> aspirin; warfarin tautomers -> same canonical form; `to_mol()` returns Mol object; None/invalid SMILES handled; works without RDKit (raises OptionalDependencyError).

---

## CF-E04: Substructure Search (SMARTS)

**Goal**: Add SMARTS-based substructure searching on individual compounds and collections.

**Files to create/modify**:
- MODIFY: `src/chemfuse/analyze/similarity.py` (add `substructure_search`, `substructure_filter`)
- MODIFY: `src/chemfuse/analyze/__init__.py` (export)
- MODIFY: `src/chemfuse/models/collection.py` (add `filter_by_substructure()`)
- MODIFY: `tests/analyze/test_similarity.py`

**API Design**:
```python
# similarity.py additions
def substructure_search(smiles_list: list[str], smarts: str) -> list[bool]
def substructure_match_atoms(smiles: str, smarts: str) -> list[tuple[int, ...]]

# On CompoundCollection:
collection.filter_by_substructure(smarts: str) -> CompoundCollection
```

**Tests**: sulfonamide SMARTS matches sulfamethoxazole; benzene ring matches aspirin; pyridine SMARTS matches nicotine; invalid SMARTS raises ValueError; collection filter returns correct subset; empty collection handled.

---

## CF-E05: Synthetic Accessibility Score (SA_Score + Fsp3)

**Goal**: Add synthetic accessibility scoring and fraction sp3 carbon metric to the compound profiling pipeline.

**Files to create/modify**:
- NEW: `src/chemfuse/compute/synthetic.py`
- MODIFY: `src/chemfuse/compute/__init__.py` (export)
- MODIFY: `src/chemfuse/models/compound.py` (add `sa_score`, `fsp3` to CompoundProperties)
- NEW: `tests/compute/test_synthetic.py`

**API Design**:
```python
# synthetic.py
def synthetic_accessibility(smiles: str) -> float | None   # 1.0 (easy) to 10.0 (hard)
def fraction_sp3(smiles: str) -> float | None              # 0.0 to 1.0
def np_likeness(smiles: str) -> float | None               # natural product-likeness score
```

**RDKit modules**: `rdkit.Chem.RDConfig` + `SA_Score.sascorer`, `rdMolDescriptors.CalcFractionCSP3`

**Tests**: ethanol SA < 3 (easy); complex natural product SA > 5; aspirin Fsp3 ~ 0.1 (mostly aromatic); cyclohexane Fsp3 = 1.0; invalid SMILES returns None; OptionalDependencyError without RDKit.

---

## CF-E06: MaxMin Diversity Picking

**Goal**: Add diversity-based compound selection for library design workflows.

**Files to create/modify**:
- NEW: `src/chemfuse/analyze/diversity.py`
- MODIFY: `src/chemfuse/analyze/__init__.py` (export)
- MODIFY: `src/chemfuse/models/collection.py` (add `pick_diverse()` method)
- NEW: `tests/analyze/test_diversity.py`

**API Design**:
```python
# diversity.py
def maxmin_pick(smiles_list: list[str], n_pick: int, fp_type: str = "morgan") -> list[int]
def diversity_score(smiles_list: list[str], fp_type: str = "morgan") -> float  # mean pairwise distance

# On CompoundCollection:
collection.pick_diverse(n: int, fp_type: str = "morgan") -> CompoundCollection
collection.diversity_score(fp_type: str = "morgan") -> float
```

**RDKit modules**: `rdkit.SimDivFilters.rdSimDivPickers.MaxMinPicker`, `DataStructs.BulkTanimotoSimilarity`

**Tests**: pick 5 from 20 returns 5 distinct indices; picks are maximally diverse (pairwise distance > random subset); pick n > len returns all; diversity_score for identical compounds = 0; collection method returns correct CompoundCollection.

---

## CF-E07: R-Group Decomposition + SAR Table

**Goal**: Given a core scaffold (SMARTS), decompose a congeneric series into R-groups and generate a structured SAR table with activity values.

**Files to create/modify**:
- NEW: `src/chemfuse/analyze/rgroup.py`
- MODIFY: `src/chemfuse/analyze/__init__.py` (export)
- MODIFY: `src/chemfuse/models/collection.py` (add `decompose_rgroups()` method)
- NEW: `tests/analyze/test_rgroup.py`

**API Design**:
```python
# rgroup.py
def decompose_rgroups(smiles_list: list[str], core_smarts: str) -> pd.DataFrame
# Returns DataFrame with columns: smiles, Core, R1, R2, ..., plus any activity values if provided

def rgroup_sar_table(compounds: list[Compound], core_smarts: str, activity_type: str = "IC50") -> pd.DataFrame
# Returns DataFrame with: smiles, R1, R2, ..., activity_value, activity_units

# On CompoundCollection:
collection.decompose_rgroups(core_smarts: str) -> pd.DataFrame
collection.sar_table(core_smarts: str, activity_type: str = "IC50") -> pd.DataFrame
```

**RDKit modules**: `rdkit.Chem.rdRGroupDecomposition.RGroupDecomposition`

**Tests**: simple benzene core with 3 substituted benzenes -> correct R-groups; unmatched compound excluded; activity values from bioactivities joined correctly; empty list handled; invalid core SMARTS raises ValueError.

---

## CF-E08: ChEMBL Data Enrichment (max_phase, assay confidence, indications, UniChem IDs)

**Goal**: Capture critical ChEMBL metadata that is already in the API response but currently discarded, plus surface UniChem cross-reference IDs on the Compound model.

**Files to create/modify**:
- MODIFY: `src/chemfuse/models/compound.py` (add `max_phase`, `molecule_type`, `drugbank_id`, `kegg_id`, `chebi_id`)
- MODIFY: `src/chemfuse/models/bioactivity.py` (add `confidence_score`, `assay_description`, `data_validity_comment`)
- MODIFY: `src/chemfuse/sources/chembl.py` (parse max_phase, molecule_type in `_parse_molecule()`; parse confidence_score in `_parse_activity()`)
- MODIFY: `src/chemfuse/sources/unichem.py` (store all cross-ref IDs on Compound, not just resolve)
- MODIFY: `src/chemfuse/sources/opentargets.py` (extend GraphQL query with `hasBeenWithdrawn`, `blackBoxWarning`, `yearOfFirstApproval`)
- MODIFY: tests for all changed sources

**API Design**:
```python
# New Compound fields:
max_phase: int | None = None              # 0-4 from ChEMBL
molecule_type: str | None = None          # "Small molecule", "Antibody", etc.
drugbank_id: str | None = None            # from UniChem
kegg_id: str | None = None                # from UniChem
chebi_id: str | None = None               # from UniChem
is_withdrawn: bool | None = None          # from Open Targets
black_box_warning: bool | None = None     # from Open Targets

# New Bioactivity fields:
confidence_score: int | None = None       # 0-9 from ChEMBL
assay_description: str | None = None
data_validity_comment: str | None = None
```

**Tests**: ChEMBL response with max_phase=4 parsed correctly; confidence_score=9 stored; UniChem response stores drugbank_id; Open Targets query returns withdrawn flag; None values handled.

---

## CF-E09: SDF Import + Enriched SDF Export

**Goal**: Read compound sets from SDF files and export enriched SDF with all computed properties as SD tags.

**Files to create/modify**:
- NEW: `src/chemfuse/core/sdf.py` (import/export)
- MODIFY: `src/chemfuse/core/export.py` (enrich export_sdf with all properties)
- MODIFY: `src/chemfuse/core/__init__.py` (export)
- MODIFY: `src/chemfuse/core/batch.py` (support SDF input in batch_search)
- NEW: `tests/core/test_sdf.py`

**API Design**:
```python
# sdf.py
def read_sdf(path: str) -> list[Compound]
def write_sdf(compounds: list[Compound], path: str, include_properties: bool = True) -> str

# batch.py addition:
# batch_search now accepts .sdf files in addition to .csv/.txt

# export.py enhancement:
# export_sdf now writes all compound properties, descriptors, druglikeness, and identifiers as SD tags
```

**RDKit modules**: `rdkit.Chem.SDMolSupplier`, `rdkit.Chem.SDWriter`, `mol.SetProp()`

**Tests**: read SDF with 3 compounds returns 3 Compound objects; SMILES extracted correctly; SD properties mapped; write SDF includes all properties as tags; round-trip (read->write->read) preserves data; empty SDF handled; file not found raises appropriate error.

---

## CF-E10: SAR Landscape Visualization

**Goal**: Create publication-quality SAR landscape plots combining chemical space coordinates, activity coloring, and activity cliff edge highlighting.

**Files to create/modify**:
- NEW: `src/chemfuse/analyze/landscape.py`
- MODIFY: `src/chemfuse/analyze/__init__.py` (export)
- MODIFY: `src/chemfuse/models/collection.py` (add `plot_sar_landscape()` method)
- NEW: `tests/analyze/test_landscape.py`

**API Design**:
```python
# landscape.py
def sar_landscape(
    smiles_list: list[str],
    activities: list[float],
    method: str = "umap",              # umap, tsne, pca
    fp_type: str = "morgan",
    cliff_threshold: float = 0.8,       # Tanimoto threshold for cliff edges
    activity_diff_threshold: float = 1.0, # minimum pIC50 difference for cliff
    colormap: str = "RdYlGn_r",
) -> dict  # returns {coords: np.ndarray, cliff_pairs: list, figure: plotly.Figure | None}

# On CompoundCollection:
collection.plot_sar_landscape(
    activity_col: str = "pic50",
    method: str = "umap",
    show_cliffs: bool = True,
) -> dict
```

**Dependencies**: plotly (optional, for interactive figures), matplotlib (fallback for static)

**Tests**: 20 compounds with activities -> returns coords (20x2) + cliff_pairs list; cliff_pairs only includes pairs above threshold; figure object is created when plotly available; handles missing activities gracefully; collection method integrates with bioactivity data.
