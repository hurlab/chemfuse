# Data Sources Documentation

Complete reference for all ChemFuse data sources, including URLs, rate limits, licenses, and data coverage.

---

## PubChem (pubchem)

### Overview

PubChem is a free, public chemical database. It aggregates data from multiple sources including FDA, NIH, and depositor submissions. Operated by the National Center for Biotechnology Information (NCBI).

### Connection Details

- **Base URL**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`
- **Protocol**: REST API (HTTPS)
- **Authentication**: Not required (public API)
- **Documentation**: https://pubchem.ncbi.nlm.nih.gov/docs/PUG-REST

### Rate Limits

- **Requests per second**: No official limit (good citizenship: ~5 req/sec)
- **Throttling**: 429 Too Many Requests if exceeded
- **Cache recommendations**: Minimum 24 hours

### Available Fields

**Compound Properties:**
- PubChem CID, IUPAC name, synonyms
- SMILES (canonical and isomeric)
- InChI, InChIKey
- Molecular formula, molecular weight
- Canonical SMILES, isotopic SMILES

**Physicochemical Properties:**
- Exact mass, monoisotopic mass
- LogP (XLogP3), surface tension
- Rotatable bonds, hydrogen bonded donors/acceptors
- TPSA, molar refractivity

**3D Structure:**
- 3D coordinates (SDF format)
- Conformer models
- Shape similarity data

**Bioassay Data:**
- PubChem BioAssay results (millions of assays)
- Bioactivity summaries
- Target information

### License

**Public Domain**: All PubChem data is in the public domain. No restrictions on use.

### Data Coverage

- **Total compounds**: 117+ million
- **Structures with 3D data**: 1.2+ million
- **Bioassay records**: 1 billion+
- **Coverage**: Small molecules, experimental compounds, drug candidates

### Special Features

- 3D coordinates for screening and visualization
- Extensive bioassay database
- FDA approval status for pharmaceuticals
- Patent information

### Query Examples

```python
import chemfuse as cf

# Search by name
results = cf.search("aspirin", source="pubchem")

# Search by SMILES
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O", source="pubchem")

# Get bioassay data
assays = aspirin.get_bioactivity(source="pubchem")
```

---

## ChEMBL (chembl)

### Overview

ChEMBL is a freely available chemical database of bioactive molecules. Provided by the European Bioinformatics Institute (EBI).

### Connection Details

- **Base URL**: `https://www.ebi.ac.uk/chembl/api/data`
- **Protocol**: REST API (HTTPS)
- **Authentication**: Not required (public API)
- **Documentation**: https://chembl.gitbook.io/chembl-interface-documentation/web-services

### Rate Limits

- **Requests per second**: 10 (strict)
- **Throttling**: 429 Too Many Requests if exceeded
- **Retry-After header**: Respected by ChemFuse client
- **Cache recommendations**: Minimum 48 hours

### Available Fields

**Compound Information:**
- ChEMBL ID, preferred name, synonyms
- SMILES (canonical), InChI
- Molecular weight, LogP
- Molecular formula

**Bioactivity Data:**
- Activity values (Ki, Kd, IC50, EC50, etc.)
- Assay type and description
- Target name and UniProt ID
- Mechanism of action
- Assay organism
- Confidence score
- Data validity comments

**Protein Targets:**
- Target name, species, organism
- UniProt ID, Ensembl ID
- Target type (protein, organism, tissue, etc.)

**Drug Information:**
- FDA approval status
- Clinical trial phase
- Indication/disease
- Formulation

### License

**Creative Commons Attribution-Share Alike 4.0**: Non-commercial and commercial use permitted with attribution.

### Data Coverage

- **Bioactive compounds**: 2+ million
- **Bioactivity records**: 20+ million
- **Targets**: 14,000+
- **Organisms**: 600+
- **Focus**: Bioactive molecules with known target/activity

### Special Features

- Extensive bioactivity data with confidence scoring
- Mechanism of action information
- Clinical development phase tracking
- FDA approval and indication data

### Query Examples

```python
import chemfuse as cf

# Search ChEMBL bioactives
results = cf.search("SARS-CoV-2", source="chembl", limit=50)

# Get bioactivity data
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O", source="chembl")
bioactivities = aspirin.get_bioactivity(source="chembl")

# Filter by mechanism of action
for record in bioactivities:
    if record.get("mechanism_of_action"):
        print(record["mechanism_of_action"])
```

---

## BindingDB (bindingdb)

### Overview

BindingDB is a public database of measured binding affinities for proteins and ligands. Provides binding data for structure-activity relationship (SAR) analysis.

### Connection Details

- **Base URL**: `https://www.bindingdb.org/axis2/services/BdbWebService`
- **Protocol**: SOAP/REST API (HTTPS)
- **Authentication**: Not required (public API)
- **Documentation**: https://www.bindingdb.org/axis2/services/BdbWebService

### Rate Limits

- **Requests per second**: 5
- **Download limit**: 10,000 results per query
- **Cache recommendations**: Minimum 7 days

### Available Fields

**Binding Affinity:**
- Ki (inhibition constant)
- Kd (dissociation constant)
- IC50 (inhibitory concentration)
- EC50 (effective concentration)
- Kmax, Tmax (pharmacokinetic parameters)

**Measurement Details:**
- Assay type and method
- Reference publication (PubMed ID)
- Measurement temperature and pH
- Binding notes and comments

**Protein Target:**
- Protein name, UniProt ID
- Target species (human, mouse, rat, etc.)
- Function/classification

### License

**CC0 1.0 Universal (Public Domain Dedication)**: No restrictions. Free to use and distribute.

### Data Coverage

- **Binding data records**: 3+ million
- **Unique protein-ligand pairs**: 1+ million
- **Proteins**: 7,000+
- **Focus**: Known binding affinities, SAR data

### Special Features

- High-quality binding affinity measurements
- Protein species information
- Publication-linked data
- Temperature and pH details

### Query Examples

```python
import chemfuse as cf

# Search BindingDB for binding data
results = cf.search("kinase inhibitor", source="bindingdb", limit=100)

# Get binding affinities
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O", source="bindingdb")
bindings = aspirin.get_bioactivity(source="bindingdb")

# Filter by binding type
for binding in bindings:
    if binding["activity_type"] == "Ki":
        print(f"Ki: {binding['activity_value']} nM")
```

---

## Open Targets (opentargets)

### Overview

Open Targets is a collaborative, open-source project that integrates genetics and safety data to support drug discovery. Provides disease-target associations from multiple evidence types.

### Connection Details

- **Base URL**: `https://api.platform.opentargets.org/api/v4`
- **Protocol**: GraphQL API (HTTPS)
- **Authentication**: Not required (public API)
- **Documentation**: https://platform-docs.opentargets.org/

### Rate Limits

- **Requests per minute**: 60
- **Query complexity**: Tracked (complex queries rate-limited more strictly)
- **Cache recommendations**: Minimum 24 hours
- **Offline data**: Monthly dataset downloads available

### Available Fields

**Target Information:**
- Target name, UniProt ID, Ensembl ID
- Target class and gene information
- Protein interactions
- Subcellular location

**Disease Information:**
- Disease name, EFO ID, ICD10 code
- Disease classification
- Parent/child disease relationships

**Association Data:**
- Overall association score (0-1 scale)
- Evidence types (genetic association, somatic mutation, RNA expression, etc.)
- Score by evidence type
- Publication count and recency

**Safety Data:**
- Adverse events by severity
- Drug contraindications
- Known toxicities

### License

**Creative Commons Attribution 4.0 International**: Non-commercial and commercial use permitted with attribution.

### Data Coverage

- **Genes/targets**: 40,000+
- **Diseases**: 25,000+
- **Associations**: 5+ million
- **Evidence types**: 20+
- **Focus**: Drug target validation, SAR, adverse event signals

### Special Features

- Integrated multi-evidence platform
- Drug safety/toxicity data
- Genetic association strength
- Publication-linked evidence

### Query Examples

```python
import chemfuse as cf

# Get disease-target associations
results = cf.search("Alzheimer's disease", source="opentargets")

# Get target associations for a compound
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O", source="opentargets")

# Cross-reference to find disease associations
xrefs = aspirin.cross_reference(target_db="opentargets")
```

---

## SureChEMBL (surecheml)

### Overview

SureChEMBL is a full-text searchable database of chemical structures extracted from millions of patent documents. Provided by Elsevier.

### Connection Details

- **Base URL**: `https://www.surechembl.org/api/data`
- **Protocol**: REST API (HTTPS)
- **Authentication**: Not required (public API)
- **Documentation**: https://www.surechembl.org/

### Rate Limits

- **Requests per second**: 3
- **Timeout**: 180 seconds per request
- **Result limit**: 100,000 per query
- **Cache recommendations**: Minimum 30 days (stable archive)

### Available Fields

**Patent Information:**
- Patent number (country and ID)
- Patent filing date, publication date
- Patent title and abstract
- Patent family information

**Chemical Structure:**
- SMILES, InChI
- Molecular formula
- Structure index/mapping

**Document Metadata:**
- Patent assignee/company
- Inventor names
- Claims containing the structure
- Full-text search

### License

**Patent Documents**: Subject to original patent licenses. Chemical information extracted for research purposes.

### Data Coverage

- **Patent documents**: 70+ million
- **Unique chemical structures**: 15+ million
- **Patent families**: 50+ years of data
- **Coverage**: Synthetic compounds, medicinal chemistry, materials

### Special Features

- Patent-associated chemistry
- Company/assignee information
- Full-text patent document search
- Historical chemistry data

### Query Examples

```python
import chemfuse as cf

# Search for compounds in patents
results = cf.search("aspirin", source="surecheml")

# Find patent information
aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O", source="surecheml")
patents = aspirin.cross_reference(target_db="surecheml")

# Check patent coverage
print(f"Found in {len(patents)} patent records")
```

---

## Cross-Reference Capabilities

### UniChem Cross-References

ChemFuse uses UniChem for mapping between databases:

```
PubChem ↔ ChEMBL ↔ BindingDB
   ↓         ↓         ↓
PubChem CID ← UniChem → ChEMBL ID
```

**UniChem coverage**: 40+ chemical databases with automated structure-based matching.

**Query example:**
```python
import chemfuse as cf

aspirin = cf.get("CC(=O)Oc1ccccc1C(=O)O")
xrefs = aspirin.cross_reference()
# Returns: {
#   "pubchem": "5102",
#   "chembl": "CHEMBL25",
#   "bindingdb": "50192168",
#   ...
# }
```

---

## Comparison Table

| Source | Focus | Records | Rate Limit | License | Coverage |
|--------|-------|---------|-----------|---------|----------|
| **PubChem** | Structures, Properties, Bioassays | 1.2B+ | Good | Public Domain | 117M compounds |
| **ChEMBL** | Bioactivity, Targets, Mechanisms | 20M+ | Good | CC-BY-SA 4.0 | 2M bioactive |
| **BindingDB** | Binding Affinity, SAR | 3M+ | Good | CC0 1.0 | 1M unique pairs |
| **Open Targets** | Disease-Target, Safety | 5M+ | Excellent | CC-BY 4.0 | 40K targets |
| **SureChEMBL** | Patents, Historical Chemistry | 15M+ | Fair | Patent-based | 70M patents |

---

## Best Practices

### Caching

Enable caching to respect rate limits:

```python
import chemfuse as cf

# Cache enabled by default with 24-hour TTL
cf.config.enable_cache = True
cf.config.cache_ttl = 86400  # 24 hours
```

### Rate Limiting

Respect source rate limits:

```python
import chemfuse as cf

# Default: 3 concurrent requests, respects per-source limits
cf.config.max_concurrent = 3
```

### Error Handling

Handle source-specific errors:

```python
import chemfuse as cf

try:
    results = cf.search("compound", source="chembl")
except cf.RateLimitError:
    print("Rate limit exceeded, retrying in 60s...")
except cf.SourceError:
    print("Source connection error")
```

### Multi-Source Enrichment

Combine data from multiple sources:

```python
import chemfuse as cf

# Search in PubChem (largest database)
results = cf.search("aspirin", source="pubchem", limit=20)

# Enrich with ChEMBL bioactivity
for compound in results:
    chembl_data = compound.cross_reference(target_db="chembl")
    if chembl_data:
        bioactivity = compound.get_bioactivity(source="chembl")
        print(f"{compound.name}: {len(bioactivity)} bioactivities")
```

---

## FAQ

**Q: Which source has the most compounds?**
A: PubChem with 117+ million compounds. ChEMBL focuses on bioactivity (2M compounds).

**Q: Can I use these sources commercially?**
A: Yes, all are CC-licensed or public domain. Check individual licenses for attribution requirements.

**Q: What's the recommended cache TTL?**
A: 24 hours for PubChem/ChEMBL, 7+ days for BindingDB, 30+ days for SureChEMBL.

**Q: Can I download bulk data?**
A: Yes, all sources provide bulk downloads. Check individual source documentation.

**Q: Which source is best for SAR analysis?**
A: BindingDB has the highest quality binding affinity data for SAR studies.
