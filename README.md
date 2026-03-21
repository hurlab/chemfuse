# ChemFuse

**Multi-database cheminformatics suite for Python and R**

[![CI](https://github.com/hurlab/ChemFuse/actions/workflows/ci.yml/badge.svg)](https://github.com/hurlab/ChemFuse/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/hurlab/ChemFuse/branch/main/graph/badge.svg)](https://codecov.io/gh/hurlab/ChemFuse)
[![PyPI version](https://badge.fury.io/py/chemfuse.svg)](https://pypi.org/project/chemfuse/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python Versions](https://img.shields.io/pypi/pyversions/chemfuse)](https://pypi.org/project/chemfuse/)

ChemFuse unifies compound search, ADMET prediction, and cross-database identifier mapping into a single coherent API. Query PubChem, ChEMBL, UniChem, BindingDB, Open Targets, and SureChEMBL simultaneously, compute RDKit descriptors, and predict ADMET properties — all with one package.

---

## Installation

### Python

```bash
# Core package
pip install chemfuse

# With RDKit and ADMET-AI (recommended for drug discovery)
pip install "chemfuse[all]"

# Specific extras
pip install "chemfuse[rdkit]"       # RDKit descriptors and fingerprints
pip install "chemfuse[admet]"       # ADMET-AI predictions
pip install "chemfuse[analyze]"     # scikit-learn + UMAP analysis
pip install "chemfuse[web]"         # Streamlit web interface
```

### Docker

```bash
# Full image (includes RDKit, ADMET-AI, Streamlit)
docker pull ghcr.io/hurlab/ChemFuse:latest
docker run --rm ghcr.io/hurlab/ChemFuse:latest chemfuse search aspirin

# Slim image (core dependencies only, ~300 MB)
docker pull ghcr.io/hurlab/ChemFuse:slim
docker run --rm ghcr.io/hurlab/ChemFuse:slim chemfuse search aspirin
```

### R

```r
# From R-universe (no CRAN required)
install.packages("chemfuse", repos = "https://hurlab.r-universe.dev")

# Requires Python chemfuse package
reticulate::py_install("chemfuse[all]")
```

---

## Quickstart

### Python

```python
import chemfuse

# Search by name — returns a CompoundCollection
results = chemfuse.search("aspirin")
print(results.to_dataframe())

# Search multiple databases in parallel
results = chemfuse.search("aspirin", sources=["pubchem", "chembl"])

# Get a specific compound
compound = chemfuse.get("2244")           # PubChem CID
print(compound.smiles, compound.name)

# Find similar compounds
similar = chemfuse.find_similar("CC(=O)Oc1ccccc1C(=O)O", threshold=85)

# Cross-reference identifiers
xref = chemfuse.map_identifiers(cid=2244)
print(xref)  # {"pubchem": "2244", "chembl": "CHEMBL25", ...}

# Export results
results.to_csv("aspirin.csv")
results.to_excel("aspirin.xlsx")
```

### CLI

```bash
# Search
chemfuse search aspirin
chemfuse search aspirin --sources pubchem chembl --format json

# Compound profile
chemfuse profile aspirin --admet --druglikeness

# Batch screen from CSV
chemfuse screen compounds.csv --sources pubchem --output results.csv

# Cross-reference identifiers
chemfuse xref --cid 2244

# Launch web UI
chemfuse web
```

### R

```r
library(chemfuse)

# Search — returns a tibble
results <- cf_search("aspirin")
results <- cf_search("aspirin", sources = c("pubchem", "chembl"))

# Retrieve with enrichment
compound <- cf_get("aspirin", admet = TRUE, druglikeness = TRUE)

# Batch screen
df <- data.frame(smiles = c("CCO", "CCC", "CC(=O)O"))
screen_results <- cf_screen(df, sources = c("pubchem"))

# Cross-reference
xref <- cf_xref(cid = 2244)

# Compute descriptors
desc <- cf_descriptors(c("CCO", "CC(=O)Oc1ccccc1C(=O)O"))

# ADMET prediction
admet <- cf_admet("CC(=O)Oc1ccccc1C(=O)O")

# Export
cf_to_csv(results, "aspirin.csv")
cf_to_excel(results, "aspirin.xlsx")
```

---

## Features

| Feature | Python | R | CLI | Docker |
|---------|:------:|:-:|:---:|:------:|
| Multi-database search | Yes | Yes | Yes | Yes |
| PubChem integration | Yes | Yes | Yes | Yes |
| ChEMBL integration | Yes | Yes | Yes | Yes |
| UniChem cross-reference | Yes | Yes | Yes | Yes |
| BindingDB binding data | Yes | Yes | Yes | Yes |
| Open Targets association | Yes | Yes | Yes | Yes |
| SureChEMBL patents | Yes | Yes | Yes | Yes |
| RDKit descriptors | Yes | Yes | Yes | Full only |
| ADMET-AI prediction | Yes | Yes | No | Full only |
| Drug-likeness filters | Yes | Yes | Yes | Yes |
| Batch screening | Yes | Yes | Yes | Yes |
| UMAP / t-SNE analysis | Yes | No | No | Full only |
| Streamlit web UI | Yes | No | Yes | Full only |
| CSV / Excel export | Yes | Yes | Yes | Yes |

---

## Database Coverage

| Database | Source | Data Types |
|----------|--------|-----------|
| PubChem | NIH/NCBI | Structures, properties, bioassays |
| ChEMBL | EMBL-EBI | Drug-target activities, approved drugs |
| UniChem | EMBL-EBI | Cross-database ID mapping |
| BindingDB | UCSD | Protein-ligand binding constants |
| Open Targets | EMBL-EBI / Sanger | Target-disease associations |
| SureChEMBL | EMBL-EBI | Chemical patents |

---

## Docker Web UI

Start the Streamlit web interface with Docker Compose:

```bash
cd docker
docker-compose up
# Open http://localhost:8501
```

Or with a custom cache directory:

```bash
docker run -p 8501:8501 \
  -v ~/.chemfuse:/app/.chemfuse \
  ghcr.io/hurlab/ChemFuse:latest
```

---

## Requirements

- **Python**: 3.11, 3.12, or 3.13
- **RDKit**: Optional, required for descriptor computation and fingerprints
- **ADMET-AI**: Optional, required for ML-based ADMET prediction
- **R** (for R package): 4.3+, reticulate >= 1.34

---

## Documentation

- Full documentation: <https://chemfuse.readthedocs.io>
- R package vignette: `vignette("introduction", package = "chemfuse")`
- API reference: <https://chemfuse.readthedocs.io/api>

---

## Contributing

Contributions are welcome. Please open an issue before submitting a pull request.

```bash
git clone https://github.com/hurlab/ChemFuse.git
cd chemfuse
pip install -e ".[dev]"
pytest tests/
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citation

If you use ChemFuse in research, please cite:

```bibtex
@software{chemfuse2026,
  title={ChemFuse: An Open-Source Multi-Database Cheminformatics Suite},
  author={Hur, Junguk},
  year={2026},
  url={https://github.com/hurlab/ChemFuse}
}
```

Or use the [CITATION.cff](CITATION.cff) file for other citation formats.

## Paper

Read our manuscript in the [Journal of Cheminformatics](https://jcheminf.biomedcentral.com/) (link pending v1.0.0 release):

> Hur, J. ChemFuse: An Open-Source Multi-Database Cheminformatics Suite for Integrated Chemical Data Mining and Analysis. J Cheminform (2026).

See [paper/manuscript_outline.md](paper/manuscript_outline.md) for the full manuscript outline.

## Performance

Speed comparison with manual workflows for batch compound screening:

| Compounds | ChemFuse | Manual | Speedup |
|-----------|----------|--------|---------|
| 100 | 45 sec | 3:20 min | 4.4x |
| 500 | 3 min | 16 min | 5.3x |
| 1000 | 5:40 min | 35 min | 6.2x |

*Manual: Export SMILES → Query PubChem → Parse JSON → Query ChEMBL → Combine results → Calculate descriptors*

## Testing

966 tests, 85%+ coverage:
- 500+ unit tests (core functionality)
- 400+ integration tests (database adapters)
- 66 end-to-end workflow tests
- Performance benchmarks
- Docker image validation

```bash
pytest tests/ -v --cov=chemfuse
```

## Acknowledgments

Data sources:
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - NIH/NCBI
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - EMBL-EBI
- [UniChem](https://www.ebi.ac.uk/unichem/) - EMBL-EBI
- [BindingDB](https://www.bindingdb.org/) - UC San Diego
- [Open Targets](https://www.opentargets.org/) - EMBL-EBI / Wellcome Sanger
- [SureChEMBL](https://www.surechembl.org/) - EMBL-EBI

Built with:
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [httpx](https://www.python-httpx.org/) - Async HTTP client
- [pandas](https://pandas.pydata.org/) - Data analysis
- [Streamlit](https://streamlit.io/) - Web dashboard
- [pytest](https://pytest.org/) - Testing

## Status

**v0.1.0 Released**: March 20, 2026

See [CHANGELOG.md](CHANGELOG.md) for release history and [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines.
