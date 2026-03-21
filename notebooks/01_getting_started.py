"""
Jupyter Notebook: Getting Started with ChemFuse

This notebook covers:
1. Installation and basic imports
2. Searching for compounds
3. Examining compound properties
4. Visualizing molecular structures

Structure:
- Cell 1: Markdown - Introduction and learning objectives
- Cell 2: Code - Install and import ChemFuse
- Cell 3: Markdown - What is ChemFuse? Overview of capabilities
- Cell 4: Code - Search for a compound (aspirin)
- Cell 5: Markdown - Understanding search results
- Cell 6: Code - Get detailed compound information
- Cell 7: Markdown - Molecular properties explained
- Cell 8: Code - Calculate molecular descriptors
- Cell 9: Markdown - Drug-likeness filters
- Cell 10: Code - Check drug-likeness criteria
- Cell 11: Markdown - Summary and next steps
- Cell 12: Code - Export results
"""

# ==============================================================================
# CELL 1: MARKDOWN - TITLE AND OBJECTIVES
# ==============================================================================
"""
# Getting Started with ChemFuse

## Learning Objectives
After completing this notebook, you will be able to:
- Install and import ChemFuse
- Search for compounds by name, SMILES, or identifier
- Retrieve detailed molecular properties
- Calculate molecular descriptors
- Evaluate drug-likeness criteria
- Export results for further analysis

## Prerequisites
- Python 3.10 or higher
- pip or conda package manager
- Basic familiarity with Python and pandas

## Time Estimate
15-20 minutes
"""

# ==============================================================================
# CELL 2: CODE - INSTALLATION AND IMPORTS
# ==============================================================================
"""
!pip install chemfuse pandas matplotlib rdkit

import chemfuse as cf
import pandas as pd
from pathlib import Path

# Verify installation
print(f"ChemFuse version: {cf.__version__}")
print(f"Available sources: {cf.AVAILABLE_SOURCES}")
"""

# ==============================================================================
# CELL 3: MARKDOWN - INTRODUCTION TO CHEMFUSE
# ==============================================================================
"""
# What is ChemFuse?

ChemFuse is an open-source cheminformatics suite that integrates multiple
chemical databases and provides tools for compound screening, analysis, and
property prediction.

## Key Capabilities
- **Multi-database search**: Search across PubChem, ChEMBL, BindingDB, and more
- **Molecular analysis**: Calculate 200+ molecular descriptors
- **Drug-likeness screening**: Apply Lipinski, Veber, and other filters
- **ADMET prediction**: Predict absorption, distribution, metabolism, excretion
- **Chemical space analysis**: Cluster and visualize compound libraries
- **Cross-database enrichment**: Find the same compound across different sources

## Data Sources
- **PubChem**: 117+ million compounds, 3D coordinates
- **ChEMBL**: 2+ million bioactive compounds, mechanism of action
- **BindingDB**: Binding affinity data (Ki, Kd, IC50)
- **Open Targets**: Disease-target associations
- **SureChEMBL**: Patent chemistry data

In this notebook, we'll focus on basic search and property analysis.
"""

# ==============================================================================
# CELL 4: CODE - SEARCH FOR A COMPOUND
# ==============================================================================
"""
# Search for aspirin in PubChem
results = cf.search("aspirin", source="pubchem", limit=10)

print(f"Found {len(results)} compounds")
print("\\nFirst 3 results:")
for i, compound in enumerate(results[:3]):
    print(f"{i+1}. {compound.name} (MW: {compound.molecular_weight:.1f})")
"""

# ==============================================================================
# CELL 5: MARKDOWN - UNDERSTANDING SEARCH RESULTS
# ==============================================================================
"""
# Understanding Search Results

The search returned a `CompoundCollection` object containing multiple compounds
that match "aspirin". Each compound in the collection has:

- **name**: IUPAC or common name
- **smiles**: Simplified Molecular Input Line Entry System (canonical representation)
- **inchi**: International Chemical Identifier (unique identifier)
- **molecular_weight**: MW in g/mol
- **formula**: Molecular formula

ChemFuse uses fuzzy matching on PubChem, so you may get results that partially
match your search term. The first result is usually the best match.
"""

# ==============================================================================
# CELL 6: CODE - GET DETAILED INFORMATION
# ==============================================================================
"""
# Get the first (best) result
aspirin = results[0]

print(f"Name: {aspirin.name}")
print(f"SMILES: {aspirin.smiles}")
print(f"InChI: {aspirin.inchi}")
print(f"Molecular Weight: {aspirin.molecular_weight} g/mol")
print(f"Molecular Formula: {aspirin.formula}")
print(f"LogP: {aspirin.logp}")
print(f"H-Bond Donors: {aspirin.hbd}")
print(f"H-Bond Acceptors: {aspirin.hba}")
print(f"Rotatable Bonds: {aspirin.rotatable_bonds}")
print(f"TPSA: {aspirin.tpsa} Ų")
"""

# ==============================================================================
# CELL 7: MARKDOWN - MOLECULAR PROPERTIES EXPLAINED
# ==============================================================================
"""
# Understanding Molecular Properties

Each compound has several important properties:

| Property | Description | Importance |
|----------|-------------|-----------|
| **Molecular Weight (MW)** | Sum of atomic masses | Important for absorption and blood-brain barrier penetration |
| **LogP** | Partition coefficient (lipophilicity) | Affects absorption, distribution, toxicity |
| **H-Bond Donors (HBD)** | Groups that can donate hydrogen bonds | Important for target binding |
| **H-Bond Acceptors (HBA)** | Groups that can accept hydrogen bonds | Important for target binding |
| **Rotatable Bonds** | Bonds allowing rotation | Affects molecular flexibility |
| **TPSA** | Topological Polar Surface Area | Predicts cell membrane permeability |

These properties are used to predict drug-likeness and ADMET (absorption,
distribution, metabolism, excretion, toxicity) properties.
"""

# ==============================================================================
# CELL 8: CODE - CALCULATE DESCRIPTORS
# ==============================================================================
"""
# Calculate all molecular descriptors
all_descriptors = aspirin.get_descriptors()

# Convert to dataframe for better viewing
descriptor_df = pd.DataFrame(
    [all_descriptors],
    index=[aspirin.name]
)

print(f"Calculated {len(all_descriptors)} descriptors")
print("\\nFirst 20 descriptors:")
print(descriptor_df.iloc[:, :20])

# Get specific descriptors
specific = aspirin.get_descriptors(['MW', 'LogP', 'TPSA', 'NumHBD', 'NumHBA'])
print("\\nLipinski-related descriptors:")
for key, value in specific.items():
    print(f"  {key}: {value}")
"""

# ==============================================================================
# CELL 9: MARKDOWN - DRUG-LIKENESS FILTERS
# ==============================================================================
"""
# Drug-Likeness Filters

Drug-likeness filters are sets of rules that predict whether a molecule is
likely to be a successful oral drug. Several filters exist:

## Lipinski's Rule of Five (most common)
- Molecular Weight ≤ 500 g/mol
- LogP ≤ 5
- H-Bond Donors ≤ 5
- H-Bond Acceptors ≤ 10

## Other Filters
- **Veber**: MW ≤ 430, TPSA ≤ 140, Rotatable bonds ≤ 10
- **Ghose**: MW 160-480, LogP -0.4 to 5.6, TPSA ≥ 20
- **Egan**: 4 ≤ LogP ≤ 5.88, 20 ≤ TPSA ≤ 130
- **Muegge**: Optimized for kinase inhibitors

These filters are heuristic rules and violations don't mean a compound can't
be a good drug, but they suggest potential issues.
"""

# ==============================================================================
# CELL 10: CODE - CHECK DRUG-LIKENESS
# ==============================================================================
"""
# Check drug-likeness for aspirin
drug_like_results = aspirin.check_drug_likeness()

print("Drug-Likeness Results for Aspirin:")
print("=" * 50)
for filter_name, result in drug_like_results.items():
    status = "✓ PASS" if result else "✗ FAIL"
    print(f"{filter_name.upper()}: {status}")

# Get details
lipinski = aspirin.check_drug_likeness(['lipinski', 'veber'])
print("\\nDetailed Lipinski Criteria:")
print(f"  MW ≤ 500: {aspirin.molecular_weight <= 500}")
print(f"  LogP ≤ 5: {aspirin.logp <= 5}")
print(f"  HBD ≤ 5: {aspirin.hbd <= 5}")
print(f"  HBA ≤ 10: {aspirin.hba <= 10}")
"""

# ==============================================================================
# CELL 11: MARKDOWN - SUMMARY AND NEXT STEPS
# ==============================================================================
"""
# Summary

Congratulations! You've completed the basics of ChemFuse:

✓ Searched for compounds in a public chemical database
✓ Retrieved molecular properties
✓ Calculated molecular descriptors
✓ Evaluated drug-likeness criteria

## Key Takeaways
1. ChemFuse integrates multiple chemical databases
2. Molecular properties predict drug-likeness
3. Drug-likeness filters are heuristic rules based on historical drug data
4. Multiple filters exist for different compound classes

## Next Steps
- Explore **Notebook 2**: Batch screening and ADMET prediction
- Explore **Notebook 3**: Multi-source enrichment workflows
- Visit the [API Reference](../docs/api-reference.md) for advanced features
- Check the [Quickstart Guide](../docs/quickstart.md) for command-line usage

## Useful Resources
- ChemFuse Documentation: https://github.com/hurlab/ChemFuse
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- RDKit (chemistry toolkit): https://www.rdkit.org/
"""

# ==============================================================================
# CELL 12: CODE - EXPORT RESULTS
# ==============================================================================
"""
# Save results for further analysis

# Convert to pandas dataframe
df = pd.DataFrame({
    'name': [c.name for c in results[:5]],
    'smiles': [c.smiles for c in results[:5]],
    'molecular_weight': [c.molecular_weight for c in results[:5]],
    'logp': [c.logp for c in results[:5]],
    'hbd': [c.hbd for c in results[:5]],
    'hba': [c.hba for c in results[:5]],
})

print("Results summary:")
print(df)

# Export to CSV
output_path = Path("aspirin_results.csv")
df.to_csv(output_path, index=False)
print(f"\\nSaved to {output_path}")
"""

if __name__ == "__main__":
    print(__doc__)
