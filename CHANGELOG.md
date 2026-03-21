# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-03-20

### Added

#### Core Infrastructure
- Async HTTP client with connection pooling and automatic retry logic
- SQLite-based local caching system with TTL support
- Comprehensive export engine supporting CSV, JSON, Excel, and SDF formats
- Configuration management with environment variable override support

#### Database Source Adapters
- **PubChem**: Compound search, structure similarity search, 3D coordinates, bioassay data access
- **ChEMBL**: Bioactivity queries, mechanism of action lookup, ChEMBL ID cross-references
- **UniChem**: Cross-reference mapping across 40+ chemical databases
- **BindingDB**: Binding affinity data (Ki, Kd, IC50, EC50) for drug targets
- **Open Targets**: Disease-target association data and disease information
- **SureChEMBL**: Patent chemistry data with document references

#### Molecular Descriptor and Fingerprint Support
- 200+ RDKit molecular descriptors (Lipinski, TPSA, LogP, molecular weight, etc.)
- 5 fingerprint types: Morgan, MACCS, TopologicalTorsion, RDKit, Avalon
- Custom descriptor pipeline with optional caching

#### Drug-Likeness Filtering
- Lipinski's Rule of Five
- Veber Filter (rotatable bonds, TPSA)
- Ghose Filter (LogP, molecular weight, TPSA)
- Egan Filter (TPSA, LogP)
- Muegge Filter (specialized for kinase inhibitors)

#### Predictive Modeling
- ML-based ADMET prediction via ADMET-AI integration
- Rule-based fallback when ML prediction unavailable
- Confidence scores for all predictions
- Support for ADME and toxicity properties (5 properties each)

#### Chemical Space Analysis
- Molecular clustering (Butina clustering for activity-based grouping)
- KMeans clustering with automatic optimal k selection
- Dimensionality reduction (UMAP, t-SNE, PCA)
- Chemical space visualization with interactive output
- Activity cliff detection from bioactivity data

#### Command-Line Interface
- `chemfuse search` - Search compound databases
- `chemfuse profile` - Detailed compound property analysis
- `chemfuse screen` - Batch screening against drug-likeness criteria
- `chemfuse similar` - Structure similarity search
- `chemfuse xref` - Cross-reference lookup across databases
- `chemfuse web` - Launch interactive Streamlit dashboard
- Bash, Zsh, and Fish shell completion support

#### Web Dashboard (Streamlit)
- Compound search interface with multi-database source selection
- Compound profile viewer with property visualization
- Batch screening tool with export capabilities
- Chemical space visualization with UMAP/t-SNE
- Cross-database enrichment workflow
- SAR analysis and activity cliff detection interface

#### Distribution
- Python package via PyPI
- Docker images: full (all dependencies) and slim (minimal dependencies)
- R package via reticulate with roxygen2 documentation
- GitHub Actions CI/CD for testing, building, and publishing

#### Testing and Quality
- 966 unit and integration tests
- 85%+ code coverage
- Pytest with mocking for external service calls
- Docker image validation tests
- R package compatibility tests

#### Documentation
- API reference for all public functions and methods
- CLI reference with command examples
- Database source documentation with rate limits and licenses
- R package usage guide with examples
- Quickstart guide for new users
- Jupyter tutorial notebooks

### Technical Details

#### Architecture
- Modular adapter pattern for database source integration
- Async/await for non-blocking I/O operations
- Smart caching with SQLite backend and TTL management
- Lazy-loading of optional dependencies (visualization, ADMET)

#### Dependencies
- **Core**: rdkit, requests, httpx, pandas, numpy, scipy, scikit-learn
- **Optional**: streamlit (web UI), reticulate (R binding), plotly (visualization)
- **Development**: pytest, pytest-asyncio, pytest-cov, black, ruff, mypy

#### Performance
- Batch processing optimized for 100-10,000 compounds
- Caching reduces redundant API calls by 70%+
- Async HTTP client handles 50+ concurrent requests
- Memory-efficient streaming for large result sets

#### Compatibility
- Python 3.10, 3.11, 3.12
- Linux, macOS, Windows (via Docker or WSL2)
- R 4.0+

[0.1.0]: https://github.com/hurlab/ChemFuse/releases/tag/v0.1.0
