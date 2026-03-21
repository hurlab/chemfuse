# Quick Start Guide

Get started with ChemFuse in 5 minutes.

## Installation

Install ChemFuse using pip:

```bash
pip install chemfuse
```

For Docker users:

```bash
docker run -it chemfuse/chemfuse:latest /bin/bash
```

For R users:

```R
install.packages("chemfuse", repos = "https://hurlab.r-universe.dev")
```

## 1. Your First Search (1 minute)

Search for aspirin across multiple databases:

```python
import chemfuse as cf

# Search for aspirin
results = cf.search("aspirin", sources=["pubchem", "chembl"])
print(f"Found {len(results)} compounds")
```

## 2. Compound Profile (2 minutes)

Get detailed information about a single compound:

```python
import chemfuse as cf

# Get compound by SMILES
compound = cf.get("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin

# View basic properties
print(compound.name)
print(compound.molecular_weight)
print(compound.smiles)

# Get ADMET prediction
admet = compound.predict_admet()
print(f"Predicted LogP: {admet['LogP']:.2f}")

# Check drug-likeness
filters = compound.check_drug_likeness()
print(f"Lipinski Pass: {filters['lipinski']}")
```

## 3. Batch Screening (2 minutes)

Screen multiple compounds for drug-likeness:

```python
import chemfuse as cf

# Create a list of SMILES strings
smiles_list = [
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "CC(=O)Nc1ccc(O)cc1",  # Acetaminophen
]

# Batch search
compounds = cf.batch_search(smiles_list, source="pubchem")

# Filter for drug-likeness
for compound in compounds:
    filters = compound.check_drug_likeness()
    if filters["lipinski"]:
        print(f"{compound.name}: Lipinski Pass ✓")
```

## 4. Cross-Database Enrichment (optional)

Combine data from multiple sources:

```python
import chemfuse as cf

# Get compound from PubChem
compound = cf.search("caffeine", source="pubchem")[0]

# Enrich with ChEMBL bioactivity
chembl_data = compound.cross_reference(target_db="chembl")
if chembl_data:
    print(f"ChEMBL bioactivities: {len(chembl_data)} records")

# Get UniChem cross-references
xrefs = compound.find_cross_references()
for db, identifier in xrefs.items():
    print(f"{db}: {identifier}")
```

## 5. Export Results

Save your results:

```python
import chemfuse as cf

# Search and collect results
compounds = cf.search("aspirin OR ibuprofen", source="pubchem", limit=20)

# Export to CSV
cf.export(compounds, "results.csv")

# Export to Excel with formulas
cf.export(compounds, "results.xlsx", include_descriptors=True)

# Export to JSON
cf.export(compounds, "results.json")
```

## Next Steps

- Read the [API Reference](api-reference.md) for detailed function documentation
- Explore the [CLI Reference](cli-reference.md) for command-line usage
- Check out the [Jupyter Notebooks](../notebooks/) for advanced examples
- Visit the [Web Dashboard](quickstart.md#web-dashboard) for interactive exploration

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

## Troubleshooting

**Import Error**: Make sure ChemFuse is installed: `pip install chemfuse`

**Connection Timeout**: ChemFuse uses cached data by default. Check your internet connection for first queries.

**ADMET Prediction Unavailable**: The ML model requires an internet connection. Falls back to rule-based prediction automatically.

## Support

- Report issues on [GitHub Issues](https://github.com/hurlab/ChemFuse/issues)
- Ask questions on [GitHub Discussions](https://github.com/hurlab/ChemFuse/discussions)
- Read the full [Documentation](../README.md)
