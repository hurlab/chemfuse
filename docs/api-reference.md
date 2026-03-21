# API Reference

Complete Python API reference for ChemFuse.

## Core Module (`chemfuse`)

### Search Functions

#### `search(query, sources=None, limit=10)`

Search for compounds across one or more databases.

**Parameters:**
- `query` (str): Search term (compound name, SMILES, InChI, CAS number, or database ID)
- `sources` (list[str], optional): List of source databases (`"pubchem"`, `"chembl"`, `"bindingdb"`, `"opentargets"`, `"surecheml"`). If None, searches all available sources
- `limit` (int, default=10): Maximum number of results per source

**Returns:**
- `CompoundCollection`: Collection of `Compound` objects

**Example:**
```python
import chemfuse as cf
results = cf.search("aspirin", sources=["pubchem", "chembl"], limit=20)
```

#### `batch_search(compounds, source="pubchem", limit=1)`

Search for multiple compounds efficiently.

**Parameters:**
- `compounds` (list[str]): List of compound identifiers (SMILES, names, or IDs)
- `source` (str, default="pubchem"): Primary source database
- `limit` (int, default=1): Results per compound

**Returns:**
- `CompoundCollection`: Collection of matched compounds

**Example:**
```python
smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"]
results = cf.batch_search(smiles_list)
```

#### `get(identifier, source=None)`

Retrieve a single compound by identifier.

**Parameters:**
- `identifier` (str): SMILES, InChI, CAS, or database ID
- `source` (str, optional): Specific source database. Auto-detected if None

**Returns:**
- `Compound`: Single compound object

**Example:**
```python
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O")
```

### Similarity Functions

#### `find_similar(compound, similarity_threshold=0.7, source="pubchem", limit=10)`

Find structurally similar compounds.

**Parameters:**
- `compound` (str or Compound): Reference compound (SMILES or Compound object)
- `similarity_threshold` (float, 0-1): Tanimoto similarity cutoff
- `source` (str): Source database
- `limit` (int): Maximum results

**Returns:**
- `CompoundCollection`: Similar compounds

**Example:**
```python
similar = cf.find_similar("CC(=O)Oc1ccccc1C(=O)O", similarity_threshold=0.8)
```

### Cross-Reference Functions

#### `cross_reference(compound, target_db=None)`

Find the same compound in other databases.

**Parameters:**
- `compound` (str or Compound): Reference compound
- `target_db` (str, optional): Specific target database. If None, searches all

**Returns:**
- `dict`: Mapping of database names to identifiers

**Example:**
```python
xrefs = cf.cross_reference("CC(=O)Oc1ccccc1C(=O)O")
# {"chembl": "CHEMBL25", "bindingdb": "50000001", ...}
```

### Collection Functions

#### `export(compounds, filepath, format=None, include_descriptors=False)`

Export compounds to file.

**Parameters:**
- `compounds` (CompoundCollection or list[Compound]): Compounds to export
- `filepath` (str): Output file path (extension determines format)
- `format` (str, optional): Format override (`"csv"`, `"json"`, `"excel"`, `"sdf"`)
- `include_descriptors` (bool): Include calculated descriptors

**Returns:**
- None (writes to file)

**Example:**
```python
cf.export(results, "compounds.csv", include_descriptors=True)
```

---

## Compound Class

### Properties

#### Basic Information
- `name` (str): IUPAC name
- `smiles` (str): SMILES string
- `inchi` (str): InChI string
- `molecular_weight` (float): MW in g/mol
- `formula` (str): Molecular formula

#### Physicochemical Properties
- `logp` (float): Partition coefficient
- `hbd` (int): Hydrogen bond donors
- `hba` (int): Hydrogen bond acceptors
- `rotatable_bonds` (int): Number of rotatable bonds
- `tpsa` (float): Topological polar surface area
- `aromatic_rings` (int): Number of aromatic rings

#### Database IDs
- `pubchem_cid` (str): PubChem CID
- `chembl_id` (str): ChEMBL ID
- `bindingdb_id` (str): BindingDB ID
- `cas_number` (str): CAS Registry Number

### Methods

#### `predict_admet()`

Predict ADMET properties using ML model or rule-based fallback.

**Returns:**
- `dict`: Predicted ADME properties with confidence scores

**Keys:**
- `"Absorption"`: {score, confidence}
- `"Distribution"`: {score, confidence}
- `"Metabolism"`: {score, confidence}
- `"Excretion"`: {score, confidence}
- `"Toxicity"`: {score, confidence}

**Example:**
```python
compound = cf.get("CC(=O)Oc1ccccc1C(=O)O")
admet = compound.predict_admet()
print(admet["Absorption"]["score"])
```

#### `check_drug_likeness(filters=None)`

Evaluate compound against drug-likeness filters.

**Parameters:**
- `filters` (list[str], optional): Specific filters to check. Default is all.

**Returns:**
- `dict`: Filter results (pass/fail) with supporting metrics

**Available Filters:**
- `"lipinski"`: Lipinski's Rule of Five
- `"veber"`: Veber filter
- `"ghose"`: Ghose filter
- `"egan"`: Egan filter
- `"muegge"`: Muegge filter

**Example:**
```python
results = compound.check_drug_likeness(filters=["lipinski", "veber"])
if results["lipinski"]:
    print("Passes Lipinski's Rule of Five")
```

#### `get_descriptors(descriptor_names=None)`

Calculate molecular descriptors.

**Parameters:**
- `descriptor_names` (list[str], optional): Specific descriptors. Default calculates 200+ standard descriptors

**Returns:**
- `dict`: Descriptor name to value mapping

**Common Descriptors:**
- `"MW"`: Molecular weight
- `"LogP"`: Partition coefficient
- `"TPSA"`: Topological polar surface area
- `"HBD"`: H-bond donors
- `"HBA"`: H-bond acceptors
- Full list of 200+ descriptors available in property reference

**Example:**
```python
descriptors = compound.get_descriptors(["MW", "LogP", "TPSA"])
```

#### `get_fingerprint(fingerprint_type="morgan")`

Get molecular fingerprint.

**Parameters:**
- `fingerprint_type` (str): `"morgan"`, `"maccs"`, `"topological"`, `"rdkit"`, or `"avalon"`

**Returns:**
- `np.ndarray`: Binary fingerprint vector

**Example:**
```python
fp = compound.get_fingerprint("morgan")
```

#### `cross_reference(target_db=None)`

Find this compound in other databases.

**Parameters:**
- `target_db` (str, optional): Specific target database

**Returns:**
- `dict`: Database name to identifier mapping

**Example:**
```python
xrefs = compound.cross_reference()
```

#### `get_bioactivity(target_name=None, source="chembl")`

Retrieve bioactivity data for this compound.

**Parameters:**
- `target_name` (str, optional): Filter by protein target
- `source` (str): Source database (`"chembl"`, `"bindingdb"`)

**Returns:**
- `list[dict]`: Bioactivity records with assay info

**Record Keys:**
- `"target_name"`: Protein target
- `"activity_value"`: Measured activity
- `"activity_type"`: Type of measurement (Ki, Kd, IC50, etc.)
- `"assay_type"`: Assay method

**Example:**
```python
bioactivities = compound.get_bioactivity(source="chembl")
for record in bioactivities:
    print(f"{record['target_name']}: {record['activity_type']}")
```

---

## CompoundCollection Class

### Properties

#### Collection Information
- `count` (int): Number of compounds in collection
- `sources` (list[str]): Database sources included
- `total_records` (int): Total records from all sources

### Methods

#### `__iter__()`

Iterate over compounds in collection.

**Example:**
```python
for compound in results:
    print(compound.name)
```

#### `__len__()`

Get number of compounds.

**Example:**
```python
print(f"Found {len(results)} compounds")
```

#### `filter(predicate)`

Filter compounds by function.

**Parameters:**
- `predicate` (callable): Function accepting Compound, returning bool

**Returns:**
- `CompoundCollection`: Filtered collection

**Example:**
```python
drug_like = results.filter(lambda c: c.check_drug_likeness()["lipinski"])
```

#### `cluster(method="butina", n_clusters=None)`

Cluster compounds by structure.

**Parameters:**
- `method` (str): `"butina"` or `"kmeans"`
- `n_clusters` (int, optional): For kmeans, number of clusters

**Returns:**
- `list[CompoundCollection]`: List of cluster collections

**Example:**
```python
clusters = results.cluster(method="butina")
for i, cluster in enumerate(clusters):
    print(f"Cluster {i}: {len(cluster)} compounds")
```

#### `visualize(method="umap", title=None)`

Visualize chemical space.

**Parameters:**
- `method` (str): `"umap"`, `"tsne"`, or `"pca"`
- `title` (str, optional): Plot title

**Returns:**
- `plotly.graph_objects.Figure`: Interactive plot

**Example:**
```python
fig = results.visualize(method="umap")
fig.show()
```

#### `find_activity_cliffs(bioactivity_field="activity_value")`

Detect activity cliffs in the collection.

**Parameters:**
- `bioactivity_field` (str): Field to analyze for cliffs

**Returns:**
- `list[tuple]`: Pairs of compounds with large activity differences

**Example:**
```python
cliffs = results.find_activity_cliffs()
for compound1, compound2, activity_diff in cliffs:
    print(f"Cliff: {compound1.name} vs {compound2.name} ({activity_diff})")
```

#### `to_dataframe()`

Convert to pandas DataFrame.

**Returns:**
- `pd.DataFrame`: Tabular representation

**Example:**
```python
df = results.to_dataframe()
df.to_csv("results.csv")
```

---

## Exceptions

#### `ChemFuseError`

Base exception for all ChemFuse errors.

#### `SourceError`

Source database connection or query error.

#### `ValidationError`

Invalid compound identifier or parameters.

#### `NotFoundError`

Compound not found in database.

#### `RateLimitError`

Source database rate limit exceeded.

---

## Configuration

Configure ChemFuse behavior:

```python
import chemfuse as cf

# Set cache TTL (seconds)
cf.config.cache_ttl = 3600 * 24  # 24 hours

# Enable/disable HTTP caching
cf.config.enable_cache = True

# Set request timeout (seconds)
cf.config.request_timeout = 30

# Set number of concurrent connections
cf.config.max_concurrent = 10

# Enable debug logging
cf.config.debug = False
```

---

## Examples

### Complete Workflow

```python
import chemfuse as cf

# Search
results = cf.search("p38 kinase inhibitor", limit=50)

# Filter
drug_like = results.filter(
    lambda c: c.check_drug_likeness()["lipinski"]
)

# Predict ADMET
for compound in drug_like:
    admet = compound.predict_admet()
    print(f"{compound.name}: {admet['Toxicity']['score']:.2f}")

# Cluster
clusters = drug_like.cluster(method="butina")

# Visualize
fig = drug_like.visualize(method="umap")
fig.show()

# Export
cf.export(drug_like, "results.csv", include_descriptors=True)
```
