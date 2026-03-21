# R Package Guide

Complete guide for using ChemFuse from R via reticulate.

---

## Installation

### From R-Universe

The recommended installation method:

```R
install.packages("chemfuse", repos = "https://hurlab.r-universe.dev")
```

### From GitHub

For development version:

```R
remotes::install_github("hurlab/ChemFuse", subdir = "r-package")
```

### Manual Installation

Requirements:
- R 4.0 or higher
- Python 3.10+ (ChemFuse requires Python >= 3.10)
- reticulate package

Steps:

```R
# Install reticulate
install.packages("reticulate")

# Setup Python environment
reticulate::py_install("chemfuse")

# Load package
library(chemfuse)
```

### Troubleshooting Installation

**Error: "No module named 'chemfuse'"**
```R
# Create Python environment
reticulate::py_install("chemfuse", envname = "chemfuse-env")

# Configure R to use environment
Sys.setenv(RETICULATE_PYTHON = "/path/to/chemfuse-env/bin/python")
```

**Error: "Python version too old"**
```R
# Upgrade Python
reticulate::install_python("3.12")
```

**Error: "reticulate not found"**
```R
install.packages("reticulate")
```

---

## Quick Start

### Load Library

```R
library(chemfuse)
```

### Your First Search

```R
# Search for aspirin
results <- cf_search("aspirin", limit = 20)

# View results as tibble
as_tibble(results)

# Extract compound names
results$name
```

### Get Compound Details

```R
# Get compound by SMILES
aspirin <- cf_get("CC(=O)Oc1ccccc1C(=O)O")

# View all properties
print(aspirin)

# Extract specific property
aspirin$molecular_weight
```

### Predict ADMET

```R
# Get ADMET predictions
admet <- cf_predict_admet(aspirin)

# View predictions
as_tibble(admet)

# Extract specific prediction
admet$absorption
```

---

## Core Functions

### Search Functions

#### `cf_search(query, sources = NULL, limit = 10, ...)`

Search for compounds across databases.

**Arguments:**
- `query` (character): Compound name, SMILES, InChI, or CAS number
- `sources` (character vector, optional): Specific sources (`"pubchem"`, `"chembl"`, `"bindingdb"`, `"opentargets"`, `"surecheml"`). Default: all
- `limit` (integer): Results per source (default: 10)
- `...`: Additional arguments passed to Python

**Returns:**
- `CompoundCollection`: Tibble-compatible collection of compounds

**Example:**
```R
# Search all sources
results <- cf_search("caffeine", limit = 20)

# Search specific sources
results <- cf_search("aspirin", sources = c("pubchem", "chembl"))

# Convert to tibble for dplyr operations
results_df <- as_tibble(results)
results_df %>%
  filter(molecular_weight < 400)
```

#### `cf_batch_search(compounds, source = "pubchem", ...)`

Search for multiple compounds efficiently.

**Arguments:**
- `compounds` (character vector): List of identifiers (SMILES, names, CAS)
- `source` (character): Primary source database (default: "pubchem")
- `...`: Additional arguments

**Returns:**
- `CompoundCollection`: Collection of matched compounds

**Example:**
```R
# Search multiple compounds
smiles <- c(
  "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
)
results <- cf_batch_search(smiles)
```

#### `cf_get(identifier, source = NULL, ...)`

Get single compound details.

**Arguments:**
- `identifier` (character): SMILES, name, CAS, or database ID
- `source` (character, optional): Specific source. Auto-detected if NULL
- `...`: Additional arguments

**Returns:**
- `Compound`: Single compound object

**Example:**
```R
aspirin <- cf_get("CC(=O)Oc1ccccc1C(=O)O")
print(aspirin)
```

### Compound Analysis Functions

#### `cf_predict_admet(compound, ...)`

Predict ADMET properties.

**Arguments:**
- `compound` (Compound or character): Compound or SMILES string
- `...`: Additional arguments

**Returns:**
- `list`: Named list of ADME properties and confidence scores

**Example:**
```R
compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")
admet <- cf_predict_admet(compound)

# Convert to tibble
as_tibble(admet)

# Check specific prediction
cat("Absorption Score:", admet$absorption$score, "\n")
cat("Confidence:", admet$absorption$confidence, "\n")
```

#### `cf_check_drug_likeness(compound, filters = NULL, ...)`

Evaluate drug-likeness filters.

**Arguments:**
- `compound` (Compound or character): Compound or SMILES string
- `filters` (character vector, optional): Specific filters (`"lipinski"`, `"veber"`, `"ghose"`, `"egan"`, `"muegge"`). Default: all
- `...`: Additional arguments

**Returns:**
- `list`: Filter results with pass/fail and supporting metrics

**Example:**
```R
compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")
results <- cf_check_drug_likeness(compound)

# Check Lipinski pass
results$lipinski

# Check specific filters
filters <- cf_check_drug_likeness(compound,
                                   filters = c("lipinski", "veber"))
```

#### `cf_get_descriptors(compound, descriptor_names = NULL, ...)`

Calculate molecular descriptors.

**Arguments:**
- `compound` (Compound or character): Compound or SMILES string
- `descriptor_names` (character vector, optional): Specific descriptors. Default: all 200+
- `...`: Additional arguments

**Returns:**
- `tibble`: Descriptor names and values

**Example:**
```R
compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")

# Get all descriptors
descriptors <- cf_get_descriptors(compound)

# Get specific descriptors
descriptors <- cf_get_descriptors(compound,
                                   c("MW", "LogP", "TPSA", "HBD", "HBA"))

# Convert to tibble
as_tibble(descriptors)
```

#### `cf_get_fingerprint(compound, fingerprint_type = "morgan", ...)`

Get molecular fingerprint.

**Arguments:**
- `compound` (Compound or character): Compound or SMILES string
- `fingerprint_type` (character): Type (`"morgan"`, `"maccs"`, `"topological"`, `"rdkit"`, `"avalon"`). Default: "morgan"
- `...`: Additional arguments

**Returns:**
- `numeric vector`: Binary fingerprint

**Example:**
```R
compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")

# Get Morgan fingerprint
fp <- cf_get_fingerprint(compound, "morgan")

# Get MACCS fingerprint
fp_maccs <- cf_get_fingerprint(compound, "maccs")

# Use for similarity calculations
fp1 <- cf_get_fingerprint(cf_get("aspirin"))
fp2 <- cf_get_fingerprint(cf_get("ibuprofen"))
similarity <- sum(fp1 * fp2) / (sqrt(sum(fp1)) * sqrt(sum(fp2)))
```

### Screening Functions

#### `cf_screen(compounds, filters = c("lipinski"), ...)`

Screen multiple compounds.

**Arguments:**
- `compounds` (character vector or CompoundCollection): Compounds to screen
- `filters` (character vector): Filters to apply (default: "lipinski")
- `...`: Additional arguments

**Returns:**
- `CompoundCollection`: Filtered compounds

**Example:**
```R
# Search candidates
candidates <- cf_search("kinase inhibitor", limit = 100)

# Screen for drug-likeness
pass_screening <- cf_screen(candidates,
                             filters = c("lipinski", "veber"))

# View results
nrow(pass_screening)
as_tibble(pass_screening)
```

### Cross-Reference Functions

#### `cf_cross_reference(compound, target_db = NULL, ...)`

Find compound in other databases.

**Arguments:**
- `compound` (Compound or character): Compound or SMILES string
- `target_db` (character, optional): Specific target database
- `...`: Additional arguments

**Returns:**
- `list`: Database to identifier mapping

**Example:**
```R
compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")

# Get all cross-references
xrefs <- cf_cross_reference(compound)
as_tibble(xrefs)

# Get specific database reference
chembl_id <- cf_cross_reference(compound, "chembl")
```

---

## Working with Results

### Convert to Tibble

All results can be converted to tibbles for dplyr operations:

```R
library(dplyr)

results <- cf_search("aspirin", limit = 50)

# Convert to tibble
df <- as_tibble(results)

# Use dplyr operations
df %>%
  filter(molecular_weight < 400) %>%
  select(name, smiles, molecular_weight) %>%
  arrange(molecular_weight)
```

### Export Results

Export to standard formats:

```R
results <- cf_search("aspirin", limit = 50)

# Export to CSV
write.csv(as_tibble(results), "results.csv")

# Export to Excel (requires openxlsx)
openxlsx::write.xlsx(as_tibble(results), "results.xlsx")

# Export to JSON
jsonlite::write_json(as_tibble(results), "results.json")
```

### Batch Processing

Process large compound sets:

```R
library(dplyr)

# Read compounds from CSV
compounds <- read.csv("compounds.csv")

# Screen each compound
results <- compounds %>%
  rowwise() %>%
  mutate(
    properties = list(cf_get_descriptors(smiles)),
    drug_like = cf_check_drug_likeness(smiles)
  ) %>%
  ungroup()

# Filter results
pass <- results %>%
  filter(drug_like$lipinski)

# Save
write.csv(pass, "pass_screening.csv")
```

---

## Advanced Usage

### Custom Descriptors

Calculate specific descriptors:

```R
# Get Lipinski-related descriptors
lipinski_descriptors <- c("MW", "LogP", "HBD", "HBA", "RotBonds")

compound <- cf_get("CC(=O)Oc1ccccc1C(=O)O")
descriptors <- cf_get_descriptors(compound, lipinski_descriptors)

# Check Lipinski rule
mw <- descriptors$MW
logp <- descriptors$LogP
hbd <- descriptors$HBD
hba <- descriptors$HBA

passes_lipinski <- (mw <= 500) & (logp <= 5) & (hbd <= 5) & (hba <= 10)
```

### Similarity Calculations

Find similar compounds:

```R
library(proxy)

# Get fingerprints
fp1 <- cf_get_fingerprint(cf_get("aspirin"))
fp2 <- cf_get_fingerprint(cf_get("ibuprofen"))
fp3 <- cf_get_fingerprint(cf_get("acetaminophen"))

# Create distance matrix
fps <- rbind(fp1, fp2, fp3)
sim_matrix <- as.matrix(dist(fps, method = "binary"))

# Find similarity
colnames(sim_matrix) <- rownames(sim_matrix) <- c("aspirin", "ibuprofen", "acetaminophen")
print(1 - sim_matrix)  # Convert distance to similarity
```

### Pipeline Integration

Create analysis pipelines:

```R
library(magrittr)

analyze_compound <- function(smiles) {
  compound <- cf_get(smiles)

  list(
    name = compound$name,
    molecular_weight = compound$molecular_weight,
    descriptors = cf_get_descriptors(compound, c("LogP", "TPSA")),
    drug_likeness = cf_check_drug_likeness(compound),
    admet = cf_predict_admet(compound)
  )
}

# Apply to multiple compounds
smiles_list <- c(
  "CC(=O)Oc1ccccc1C(=O)O",
  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
)

results <- lapply(smiles_list, analyze_compound)
```

---

## Troubleshooting

### Import Errors

**Error: `Error in library(chemfuse): there is no package called 'chemfuse'`**

Install reticulate and Python package:
```R
install.packages("reticulate")
reticulate::py_install("chemfuse")
```

### Python Environment Issues

**Error: `ModuleNotFoundError: No module named 'chemfuse'`**

Create dedicated environment:
```R
reticulate::virtualenv_create("chemfuse-env")
reticulate::virtualenv_install("chemfuse-env", "chemfuse")
reticulate::use_virtualenv("chemfuse-env")
```

### Connection Errors

**Error: `Connection timeout` or `HTTPError`**

Set timeout and retry:
```R
# Increase timeout
options(reticulate.python_version_override = "3.12")

# Check Python setup
reticulate::py_config()
```

### Data Type Mismatches

**Error: `cannot coerce class c(...) to a data frame`**

Convert to tibble explicitly:
```R
results <- cf_search("aspirin")
df <- as_tibble(results)
```

---

## FAQ

**Q: Can I use ChemFuse without reticulate?**
A: No, the R package requires reticulate to interface with Python. This is the recommended approach.

**Q: How do I use a specific Python version?**
A: Set `RETICULATE_PYTHON` environment variable:
```R
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
```

**Q: Can I install without internet?**
A: Yes, install ChemFuse Python package first, then R package from local files.

**Q: How do I speed up batch processing?**
A: Use caching:
```R
cf_config(enable_cache = TRUE)
cf_config(cache_ttl = 86400)  # 24 hours
```

**Q: How do I debug issues?**
A: Enable Python output:
```R
reticulate::py_capture_output()
```

---

## Support

- Report issues on [GitHub Issues](https://github.com/hurlab/ChemFuse/issues)
- Ask questions on [GitHub Discussions](https://github.com/hurlab/ChemFuse/discussions)
- Check the [Python documentation](api-reference.md) for underlying API details
