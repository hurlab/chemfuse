# ChemFuse MCP Server

Exposes ChemFuse cheminformatics functionality as MCP (Model Context Protocol) tools
so AI agents (Claude, Cursor, etc.) can search chemical databases, predict ADMET
properties, evaluate drug-likeness, and analyse molecular structure in natural language
workflows.

---

## Install

```bash
pip install chemfuse[mcp]
```

This installs the `mcp` SDK alongside ChemFuse. RDKit is required for local
computation tools (ADMET, drug-likeness, descriptors, scaffold analysis):

```bash
pip install chemfuse[mcp,rdkit]
```

---

## Run the server

```bash
python -m chemfuse.mcp
```

The server communicates over stdio (MCP standard transport).

---

## Configure in Claude Desktop

Add the following to `~/Library/Application Support/Claude/claude_desktop_config.json`
(macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "chemfuse": {
      "command": "python",
      "args": ["-m", "chemfuse.mcp"]
    }
  }
}
```

If you use a virtual environment, specify the full path to the interpreter:

```json
{
  "mcpServers": {
    "chemfuse": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "chemfuse.mcp"]
    }
  }
}
```

---

## Configure in Claude Code

Add the following to `.claude/mcp.json` in your project root (or `~/.claude/mcp.json`
for global access):

```json
{
  "mcpServers": {
    "chemfuse": {
      "command": "python",
      "args": ["-m", "chemfuse.mcp"],
      "env": {}
    }
  }
}
```

---

## Available tools

### `search_compounds`

Search for chemical compounds across databases.

| Parameter    | Type          | Default      | Description                               |
|--------------|---------------|--------------|-------------------------------------------|
| `query`      | string        | **required** | Name, SMILES, CID, formula, or InChI      |
| `sources`    | list[string]  | `["pubchem"]`| Databases to query (`pubchem`, `chembl`)  |
| `query_type` | string        | `"name"`     | `name`, `smiles`, `cid`, `formula`, `inchi` |
| `limit`      | integer       | `10`         | Max results (up to 100)                   |

**Example prompt:**
> "Search for compounds similar in name to ibuprofen in PubChem"

---

### `get_compound`

Retrieve a single compound by its database identifier.

| Parameter    | Type   | Default      | Description                       |
|--------------|--------|--------------|-----------------------------------|
| `identifier` | string | **required** | Database identifier (e.g. CID)    |
| `source`     | string | `"pubchem"`  | Source database                   |

**Example prompt:**
> "Get the full details of PubChem CID 2244"

---

### `find_similar`

Find structurally similar compounds using Tanimoto fingerprint similarity.

| Parameter    | Type    | Default | Description                         |
|--------------|---------|---------|-------------------------------------|
| `smiles`     | string  | **required** | SMILES of the query compound   |
| `threshold`  | integer | `90`    | Tanimoto threshold (0-100)          |
| `max_results`| integer | `10`    | Maximum number of results           |

**Example prompt:**
> "Find compounds similar to aspirin (CC(=O)Oc1ccccc1C(=O)O) with at least 80% similarity"

---

### `cross_reference`

Map a compound identifier across multiple chemical databases via UniChem.

| Parameter         | Type   | Default      | Description                                          |
|-------------------|--------|--------------|------------------------------------------------------|
| `identifier`      | string | **required** | The identifier value                                 |
| `identifier_type` | string | **required** | One of `cid`, `chembl_id`, `inchikey`, `smiles`      |

**Example prompt:**
> "Cross-reference PubChem CID 2244 to find its ChEMBL and DrugBank IDs"

---

### `predict_admet`

Predict ADMET properties for a molecule. Uses ML models (admet-ai) when installed,
falls back to rule-based RDKit heuristics otherwise.

| Parameter | Type   | Default      | Description         |
|-----------|--------|--------------|---------------------|
| `smiles`  | string | **required** | SMILES of the molecule |

**Returns:** Predictions for solubility, GI absorption, BBB permeability, CYP1A2/2C9/2C19/2D6/3A4
inhibition, P-gp substrate, hERG liability, AMES mutagenicity, hepatotoxicity, DILI,
clearance, half-life, Caco-2 permeability, HIA, and an `overall_score` (0-1).

**Example prompt:**
> "Predict ADMET properties for aspirin: CC(=O)Oc1ccccc1C(=O)O"

---

### `check_drug_likeness`

Run all standard drug-likeness filters on a molecule.

| Parameter | Type   | Default      | Description         |
|-----------|--------|--------------|---------------------|
| `smiles`  | string | **required** | SMILES of the molecule |

**Returns:** Pass/fail and violations for Lipinski Ro5, Veber, Ghose, Egan, Muegge,
PAINS, and QED score.

**Example prompt:**
> "Does this compound pass Lipinski's Rule of Five: CC(=O)Oc1ccccc1C(=O)O"

---

### `compute_descriptors`

Compute physicochemical molecular descriptors using RDKit.

| Parameter | Type   | Default      | Description         |
|-----------|--------|--------------|---------------------|
| `smiles`  | string | **required** | SMILES of the molecule |

**Returns:** Key descriptors (MW, LogP, TPSA, HBD, HBA, rotatable bonds, ring counts,
QED, etc.) with the total count of all 200+ RDKit descriptors computed.

**Example prompt:**
> "What are the molecular descriptors for caffeine: Cn1c(=O)c2c(ncn2C)n(c1=O)C"

---

### `analyze_scaffolds`

Identify and count Bemis-Murcko scaffolds across a list of compounds.

| Parameter     | Type          | Default      | Description                  |
|---------------|---------------|--------------|------------------------------|
| `smiles_list` | list[string]  | **required** | List of SMILES to analyse    |

**Returns:** Scaffold SMILES with frequency counts and relative frequencies, sorted by
most common scaffold.

**Example prompt:**
> "Analyse the scaffolds in these 5 compounds: [list of SMILES]"

---

### `compare_compounds`

Compare two compounds on similarity, properties, and scaffold identity.

| Parameter  | Type   | Default      | Description              |
|------------|--------|--------------|--------------------------|
| `smiles_a` | string | **required** | SMILES of first compound |
| `smiles_b` | string | **required** | SMILES of second compound|

**Returns:** Tanimoto similarity (with high/medium/low label), property comparison
(MW, LogP, TPSA, HBD, HBA, rotatable bonds), and shared Murcko scaffold if any.

**Example prompt:**
> "Compare aspirin and ibuprofen: how similar are they structurally?"

---

## Notes

- Network tools (`search_compounds`, `get_compound`, `find_similar`, `cross_reference`)
  make live API calls to PubChem and other databases.
- Computation tools (`predict_admet`, `check_drug_likeness`, `compute_descriptors`,
  `analyze_scaffolds`, `compare_compounds`) run entirely locally using RDKit.
- The `mcp` package is optional: `import chemfuse.mcp.server` succeeds even without it,
  so the rest of ChemFuse remains importable in environments where `mcp` is not installed.
