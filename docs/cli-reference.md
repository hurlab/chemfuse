# CLI Reference

Complete reference for ChemFuse command-line interface.

## Installation

Shell completion is automatically available after installation. For manual setup:

```bash
# Bash
eval "$(_CHEMFUSE_COMPLETE=bash_source chemfuse)"

# Zsh
eval "$(_CHEMFUSE_COMPLETE=zsh_source chemfuse)"

# Fish
_CHEMFUSE_COMPLETE=fish_source chemfuse | source
```

---

## Commands

### `chemfuse search`

Search for compounds across multiple databases.

**Usage:**
```bash
chemfuse search QUERY [OPTIONS]
```

**Arguments:**
- `QUERY`: Search term (compound name, SMILES, InChI, CAS number, or database ID)

**Options:**
- `--source TEXT`: Database source (`pubchem`, `chembl`, `bindingdb`, `opentargets`, `surecheml`, or comma-separated list). Default: all
- `--limit INTEGER`: Maximum results per source (default: 10)
- `--output TEXT`: Output file (CSV, JSON, Excel, SDF). Default: stdout
- `--format TEXT`: Override output format (`csv`, `json`, `excel`, `sdf`)
- `--include-descriptors`: Include calculated molecular descriptors
- `--include-bioactivity`: Include bioactivity data from ChEMBL
- `-v, --verbose`: Show detailed output

**Examples:**

Search all sources:
```bash
chemfuse search aspirin --limit 20
```

Search specific sources with descriptors:
```bash
chemfuse search "CC(=O)Oc1ccccc1C(=O)O" --source pubchem,chembl \
  --include-descriptors --output results.csv
```

---

### `chemfuse profile`

Show detailed properties for a single compound.

**Usage:**
```bash
chemfuse profile IDENTIFIER [OPTIONS]
```

**Arguments:**
- `IDENTIFIER`: SMILES string, compound name, CAS number, or database ID

**Options:**
- `--source TEXT`: Specify compound source. Default: auto-detect
- `--descriptors`: Calculate all 200+ molecular descriptors
- `--fingerprints`: Generate fingerprints (morgan, maccs, topological, rdkit, avalon)
- `--bioactivity`: Retrieve bioactivity data from ChEMBL
- `--drug-likeness`: Check drug-likeness filters (Lipinski, Veber, Ghose, Egan, Muegge)
- `--admet`: Predict ADMET properties
- `--format TEXT`: Output format (`text`, `json`, `csv`). Default: text
- `-v, --verbose`: Detailed output

**Examples:**

Show basic properties:
```bash
chemfuse profile "CC(=O)Oc1ccccc1C(=O)O"
```

Complete analysis with predictions:
```bash
chemfuse profile "CC(=O)Oc1ccccc1C(=O)O" \
  --descriptors --drug-likeness --admet --bioactivity
```

Output as JSON:
```bash
chemfuse profile aspirin --format json > aspirin.json
```

---

### `chemfuse screen`

Batch screen compounds against drug-likeness criteria.

**Usage:**
```bash
chemfuse screen INPUT [OPTIONS]
```

**Arguments:**
- `INPUT`: Input file (CSV, SDF, JSON, or SMILES list)

**Options:**
- `--smiles-column TEXT`: CSV column name containing SMILES. Default: "smiles"
- `--name-column TEXT`: CSV column name for compound names. Default: "name"
- `--filters TEXT`: Filters to apply (comma-separated: `lipinski`, `veber`, `ghose`, `egan`, `muegge`). Default: lipinski
- `--predict-admet`: Predict ADMET properties for passing compounds
- `--cluster`: Cluster results by structure
- `--output TEXT`: Output file (CSV, JSON, SDF)
- `--output-format TEXT`: Override output format
- `-v, --verbose`: Show screening progress

**Examples:**

Basic screening (Lipinski filter):
```bash
chemfuse screen compounds.csv --output results.csv
```

Screen with multiple filters:
```bash
chemfuse screen compounds.csv \
  --filters lipinski,veber,egan \
  --smiles-column "SMILES" \
  --output pass.csv
```

Screen with ADMET prediction:
```bash
chemfuse screen compounds.csv \
  --predict-admet \
  --output pass_admet.csv
```

---

### `chemfuse similar`

Find structurally similar compounds.

**Usage:**
```bash
chemfuse similar COMPOUND [OPTIONS]
```

**Arguments:**
- `COMPOUND`: Reference compound (SMILES, name, or CAS number)

**Options:**
- `--threshold FLOAT`: Tanimoto similarity threshold (0.0-1.0, default: 0.7)
- `--fingerprint TEXT`: Fingerprint type (`morgan`, `maccs`, `topological`, `rdkit`, `avalon`). Default: morgan
- `--source TEXT`: Search source. Default: pubchem
- `--limit INTEGER`: Maximum results (default: 20)
- `--output TEXT`: Output file (CSV, JSON, SDF)
- `-v, --verbose`: Show detailed output

**Examples:**

Find similar compounds with default settings:
```bash
chemfuse similar "CC(=O)Oc1ccccc1C(=O)O"
```

Strict similarity search:
```bash
chemfuse similar "CC(=O)Oc1ccccc1C(=O)O" \
  --threshold 0.9 \
  --fingerprint maccs \
  --limit 50 \
  --output similar.csv
```

---

### `chemfuse xref`

Cross-reference compound across databases.

**Usage:**
```bash
chemfuse xref COMPOUND [OPTIONS]
```

**Arguments:**
- `COMPOUND`: Compound identifier (SMILES, name, CAS, or database ID)

**Options:**
- `--source TEXT`: Source database for lookup
- `--target-db TEXT`: Specific target database. Default: all
- `--output TEXT`: Output file (JSON recommended)
- `--format TEXT`: Output format (`json`, `csv`, `table`). Default: table

**Examples:**

Find all cross-references:
```bash
chemfuse xref "CC(=O)Oc1ccccc1C(=O)O"
```

Output example:
```
Compound: Aspirin

Cross-references:
  pubchem:     5102
  chembl:      CHEMBL25
  bindingdb:   50192168
  opentargets: CHEMBL25
  surecheml:   X0000123456
```

Get cross-references as JSON:
```bash
chemfuse xref aspirin --format json > xrefs.json
```

---

### `chemfuse web`

Launch interactive Streamlit web dashboard.

**Usage:**
```bash
chemfuse web [OPTIONS]
```

**Options:**
- `--port INTEGER`: Server port (default: 8501)
- `--host TEXT`: Server host (default: localhost)
- `--theme CHOICE`: UI theme (`light`, `dark`). Default: light
- `--no-cache`: Disable local caching
- `--no-browser`: Don't open browser automatically

**Examples:**

Launch with default settings:
```bash
chemfuse web
```

Custom host and port:
```bash
chemfuse web --host 0.0.0.0 --port 8080
```

Dark theme, no cache:
```bash
chemfuse web --theme dark --no-cache
```

Once running, visit `http://localhost:8501` to access:
- Compound search across all databases
- Property visualization and analysis
- Batch screening workflows
- Chemical space exploration (UMAP, t-SNE)
- SAR analysis tools
- Multi-database enrichment

---

## Global Options

Available with all commands:

- `--version`: Show ChemFuse version
- `--help`: Show command help
- `--config TEXT`: Custom config file path
- `--cache-dir TEXT`: Custom cache directory
- `--no-cache`: Disable caching
- `--timeout INTEGER`: Request timeout in seconds
- `--max-concurrent INTEGER`: Maximum concurrent API requests
- `--debug`: Enable debug logging

**Examples:**

```bash
# Show version
chemfuse --version

# Use custom config
chemfuse --config ~/.chemfuse.ini search aspirin

# Disable caching
chemfuse --no-cache search aspirin

# Debug logging
chemfuse --debug search aspirin
```

---

## Configuration File

Create `~/.chemfuse.ini` to configure defaults:

```ini
[general]
cache_dir = ~/.chemfuse/cache
cache_ttl = 86400
request_timeout = 30
max_concurrent = 10
debug = false

[sources]
pubchem_enabled = true
chembl_enabled = true
bindingdb_enabled = true
opentargets_enabled = true
surecheml_enabled = true

[output]
default_format = csv
include_descriptors = false
include_bioactivity = false

[admet]
model = ml
fallback_to_rules = true
```

---

## Exit Codes

- `0`: Success
- `1`: General error
- `2`: Invalid arguments
- `3`: Source connection error
- `4`: Rate limit exceeded
- `5`: Configuration error

---

## Examples

### Complete Workflow

Search, screen, and analyze compounds:

```bash
# Search for compounds
chemfuse search "p38 kinase inhibitor" --limit 100 \
  --source chembl --output candidates.csv

# Screen for drug-likeness
chemfuse screen candidates.csv \
  --filters lipinski,veber \
  --predict-admet \
  --output pass_screening.csv

# Analyze similar compounds
chemfuse similar "CC(=O)Oc1ccccc1C(=O)O" \
  --threshold 0.8 --output similar.csv

# Cross-reference top hit
chemfuse xref "CC(=O)Oc1ccccc1C(=O)O" --format json > xrefs.json
```

### Batch Processing

Process multiple screening rounds:

```bash
# Round 1: Broad search
chemfuse search "kinase inhibitor" --limit 500 \
  --output round1_candidates.csv

# Round 2: Screen
chemfuse screen round1_candidates.csv \
  --filters lipinski,veber,ghose \
  --predict-admet \
  --output round2_screened.csv

# Round 3: Cluster
chemfuse screen round2_screened.csv \
  --cluster \
  --output round3_clustered.csv
```

### Automated Pipelines

Save as shell script for reproducible workflows:

```bash
#!/bin/bash
# compound_workflow.sh

TARGET=$1
OUTPUT_DIR=$2

chemfuse search "$TARGET" --limit 100 --output "$OUTPUT_DIR/search.csv"
chemfuse screen "$OUTPUT_DIR/search.csv" \
  --filters lipinski,veber \
  --predict-admet \
  --output "$OUTPUT_DIR/screened.csv"
echo "Complete: $OUTPUT_DIR"
```

Run with:
```bash
bash compound_workflow.sh "SARS-CoV-2" results/
```
