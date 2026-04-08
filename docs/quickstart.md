# Quick Start Guide

Get started with ChemFuse in 5 minutes.

## Installation

Install ChemFuse with all optional dependencies:

```bash
pip install chemfuse[all]
```

Or install with only the core dependencies:

```bash
pip install chemfuse
```

## 1. Your First Search (1 minute)

Search for aspirin across multiple databases:

```python
import chemfuse as cf

# Search for aspirin
results = cf.search("aspirin", sources=["pubchem", "chembl"])
print(f"Found {len(results)} compounds")

# Access the first compound
compound = results[0]
print(compound.name)
print(compound.properties.molecular_weight)
print(compound.smiles)
```

## 2. Compound Profile (2 minutes)

Get detailed information about a specific compound:

```python
import chemfuse as cf

# Get compound by PubChem CID
compound = cf.get("2244")  # Aspirin CID

# Compute molecular descriptors (requires RDKit)
descriptors = compound.compute_descriptors()
print(f"Computed {len(descriptors)} descriptors")

# Check drug-likeness (5 standard filters + PAINS + QED)
dl = compound.check_drug_likeness()
print(f"Lipinski: {'Pass' if dl.lipinski.pass_filter else 'Fail'}")
print(f"Veber: {'Pass' if dl.veber.pass_filter else 'Fail'}")
if dl.qed:
    print(f"QED Score: {dl.qed['qed']:.2f}")

# Predict ADMET properties
from chemfuse.compute.admet import predict_admet
profile = predict_admet(compound.smiles)
print(f"ADMET Score: {profile.overall_score:.2f}")
print(f"Predictions: {len(profile.predictions)}")
```

## 3. Batch Screening (2 minutes)

Screen a collection of compounds:

```python
import chemfuse as cf

# Search for multiple compounds
results = cf.search("aspirin", sources=["pubchem"])

# Compute properties for all compounds
results.compute_all(descriptors=True, druglikeness=True)

# Filter by drug-likeness
druglike = results.filter_by_druglikeness(lipinski=True, veber=True)
print(f"Drug-like compounds: {len(druglike)} / {len(results)}")

# Filter by property ranges
filtered = results.filter(mw_range=(150, 500), logp_range=(-1, 5))
print(f"Filtered: {len(filtered)} compounds")
```

## 4. Cross-Database Enrichment

Combine data from multiple sources:

```python
import chemfuse as cf

# Search and get a compound
results = cf.search("caffeine", sources=["pubchem"])
compound = results[0]

# Cross-reference across databases
enriched = cf.cross_reference(compound, target_sources=["chembl"])
print(f"ChEMBL ID: {enriched.chembl_id}")

# Map identifiers via UniChem
xrefs = cf.map_identifiers(cid=compound.cid)
for db, identifier in xrefs.items():
    print(f"{db}: {identifier}")
```

## 5. Export Results

Save your results in multiple formats:

```python
import chemfuse as cf

# Search and collect results
results = cf.search("aspirin", sources=["pubchem"])

# Export to CSV
results.to_csv("results.csv")

# Export to Excel
results.to_excel("results.xlsx")

# Export to JSON
results.to_json("results.json")

# Export to pandas DataFrame
df = results.to_dataframe()
print(df.head())
```

## Web Dashboard

Launch the interactive Streamlit dashboard:

```bash
chemfuse web
```

Then open your browser to `http://localhost:8501` to:
- Search across all databases
- Visualize compound properties
- Run batch screenings
- Explore chemical space
- Download results

## Next Steps

- Read the [API Reference](api-reference.md) for detailed function documentation
- Explore the [CLI Reference](cli-reference.md) for command-line usage
- Check out the [Jupyter Notebooks](../notebooks/) for advanced examples

## Troubleshooting

**Import Error**: Make sure ChemFuse is installed: `pip install chemfuse[all]`

**RDKit Not Found**: Install with: `pip install chemfuse[rdkit]`

**Connection Timeout**: ChemFuse caches API responses for 7 days. Check your internet connection for first queries.

**ADMET Prediction**: Uses admet-ai ML models when installed, falls back to rule-based heuristics automatically.

## Support

- Report issues on [GitHub Issues](https://github.com/hurlab/ChemFuse/issues)
- Read the full [Documentation](../README.md)
