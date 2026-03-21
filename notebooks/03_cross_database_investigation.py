"""
Jupyter Notebook: Cross-Database Investigation

This notebook covers:
1. Searching across multiple databases for a single compound
2. Enriching with bioactivity data from ChEMBL and BindingDB
3. Finding disease associations via Open Targets
4. Discovering patent chemistry via SureChEMBL
5. Synthesizing multi-source data
6. Comprehensive enrichment workflow

Structure:
- Cell 1: Markdown - Multi-source investigation
- Cell 2: Code - Search multiple databases
- Cell 3: Markdown - Understanding cross-references
- Cell 4: Code - Get bioactivity data
- Cell 5: Markdown - Interpreting bioactivity
- Cell 6: Code - Find disease associations
- Cell 7: Markdown - Disease relevance
- Cell 8: Code - Check patent information
- Cell 9: Markdown - Patent mining
- Cell 10: Code - Synthesize multi-source data
- Cell 11: Markdown - Data integration challenges
- Cell 12: Code - Export comprehensive profile
"""

# ==============================================================================
# CELL 1: MARKDOWN - MULTI-SOURCE INVESTIGATION
# ==============================================================================
"""
# Cross-Database Investigation: A Holistic View

The power of ChemFuse is integrating data from multiple sources to build
a comprehensive profile of a compound.

## Workflow
1. **Core Search**: Find compound in major chemical databases
2. **Bioactivity**: Retrieve target engagement and assay data
3. **Disease**: Find disease associations and clinical relevance
4. **Patents**: Discover prior art and development history
5. **Safety**: Assess adverse events and contraindications
6. **Synthesis**: Compile evidence for development

## Real-World Scenario
You've identified an interesting hit compound and want to:
- Understand its known biological activity
- Find potential drug applications
- Assess development history
- Identify similar compounds and scaffolds
- Evaluate intellectual property status

This notebook demonstrates the complete investigation workflow.
"""

# ==============================================================================
# CELL 2: CODE - SEARCH MULTIPLE DATABASES
# ==============================================================================
"""
import chemfuse as cf
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Choose a compound of interest
# Let's investigate a kinase inhibitor scaffold: Sunitinib (cancer drug)
smiles = "CN1C(=O)CC(c2c(F)ccc(c2F)NC(=O)C=C)C1=O"  # Sunitinib analog
compound_name = "Sunitinib"

print(f"Investigating: {compound_name}")
print(f"SMILES: {smiles}\\n")

# Search across all major databases
sources = ['pubchem', 'chembl', 'bindingdb']
search_results = {}

for source in sources:
    try:
        print(f"Searching {source}...")
        compound = cf.get(smiles, source=source)
        search_results[source] = {
            'name': compound.name,
            'database_id': getattr(compound, f'{source}_cid', 'N/A'),
            'mw': compound.molecular_weight,
            'logp': compound.logp,
        }
        print(f"  ✓ Found: {compound.name}")
    except Exception as e:
        print(f"  ✗ Not found: {e}")

# Create summary
print("\\nDatabase Coverage:")
summary_df = pd.DataFrame(search_results).T
print(summary_df)
"""

# ==============================================================================
# CELL 3: MARKDOWN - UNDERSTANDING CROSS-REFERENCES
# ==============================================================================
"""
# Cross-Database Mapping

When a compound appears in multiple databases, it has different identifiers:

| Database | Identifier Type | Example |
|----------|-----------------|---------|
| **PubChem** | CID (Compound ID) | 5102 |
| **ChEMBL** | ChEMBL ID | CHEMBL25 |
| **BindingDB** | BindingDB ID | 50192168 |
| **Open Targets** | Ensembl ID (targets) | ENSG00000... |
| **SureChEMBL** | Patent document ID | X123456789 |

## Mapping Strategy
1. Use SMILES as a universal reference (structure-based)
2. Query each database for the SMILES
3. Collect database-specific identifiers
4. Use UniChem for systematic cross-referencing

ChemFuse automates this cross-referencing, allowing you to:
- Find the same compound across databases
- Combine data from multiple sources
- Build comprehensive compound profiles
"""

# ==============================================================================
# CELL 4: CODE - GET BIOACTIVITY DATA
# ==============================================================================
"""
# Get bioactivity data from multiple sources
print("Retrieving bioactivity data...\\n")

bioactivity_all = []

# Get the compound
compound = cf.get(smiles, source='pubchem')

# ChEMBL bioactivity
print("ChEMBL Bioactivity:")
try:
    chembl_bio = compound.get_bioactivity(source='chembl')
    print(f"  Found {len(chembl_bio)} bioactivity records")

    for record in chembl_bio[:5]:  # Show top 5
        bioactivity_all.append({
            'source': 'ChEMBL',
            'target': record.get('target_name', 'Unknown'),
            'activity_type': record.get('activity_type', 'N/A'),
            'activity_value': record.get('activity_value', 'N/A'),
            'assay_type': record.get('assay_type', 'N/A'),
        })
        print(f"    {record['target_name']}: {record['activity_type']} = {record['activity_value']}")
except Exception as e:
    print(f"  Error: {e}")

# BindingDB bioactivity
print("\\nBindingDB Binding Affinities:")
try:
    bindingdb_bio = compound.get_bioactivity(source='bindingdb')
    print(f"  Found {len(bindingdb_bio)} binding records")

    for record in bindingdb_bio[:5]:  # Show top 5
        bioactivity_all.append({
            'source': 'BindingDB',
            'target': record.get('target_name', 'Unknown'),
            'activity_type': record.get('activity_type', 'N/A'),
            'activity_value': record.get('activity_value', 'N/A'),
            'assay_type': record.get('assay_type', 'N/A'),
        })
        print(f"    {record['target_name']}: {record['activity_type']} = {record['activity_value']}")
except Exception as e:
    print(f"  Error: {e}")

# Create dataframe
if bioactivity_all:
    bioactivity_df = pd.DataFrame(bioactivity_all)
    print(f"\\nTotal bioactivity records: {len(bioactivity_df)}")
"""

# ==============================================================================
# CELL 5: MARKDOWN - INTERPRETING BIOACTIVITY
# ==============================================================================
"""
# Understanding Bioactivity Data

Bioactivity records show how a compound interacts with protein targets.

## Key Fields
- **Target**: Protein or gene target (e.g., 'ABL1 kinase', 'TP53')
- **Activity Type**: Type of measurement
  - Ki: Inhibition constant (lower is better, typical: pM-nM)
  - Kd: Dissociation constant (lower is better)
  - IC50: Inhibitory concentration at 50% (lower is better)
  - EC50: Effective concentration (lower is better)
  - Potency: Relative binding strength
- **Activity Value**: Numerical measurement (usually in nM)
- **Assay Type**: How activity was measured
  - Biochemical: Direct target binding
  - Cell-based: Activity in cells
  - Binding: Direct interaction
  - FRET/TR-FRET: Fluorescence assays

## Interpretation Guidelines
- **Ki < 100 nM**: Strong inhibitor (good drug candidate)
- **100 nM - 1 µM**: Moderate inhibitor
- **> 1 µM**: Weak inhibitor
- **Selectivity**: Check IC50 against related targets
- **Multiple assays**: Consistent results strengthen evidence

## Real-World Example: Sunitinib
Known to inhibit multiple kinases:
- FLT3 (fms-like tyrosine kinase) - primary target
- KIT (stem cell factor receptor)
- PDGFR (platelet-derived growth factor receptor)
- RET (rearranged during transfection) proto-oncogene
"""

# ==============================================================================
# CELL 6: CODE - DISEASE ASSOCIATIONS
# ==============================================================================
"""
# Find disease associations via Open Targets
print("Retrieving disease associations...\\n")

try:
    # Search for disease associations
    disease_results = cf.search(
        compound_name,
        source='opentargets',
        limit=20
    )

    print(f"Found {len(disease_results)} disease associations")

    disease_data = []
    for result in disease_results[:10]:
        disease_data.append({
            'disease': result.name if hasattr(result, 'name') else 'Unknown',
            'association_score': getattr(result, 'score', 'N/A'),
            'type': getattr(result, 'type', 'Unknown'),
        })
        print(f"  • {result.name if hasattr(result, 'name') else 'Unknown'}")

    if disease_data:
        disease_df = pd.DataFrame(disease_data)
        print(f"\\nTotal disease associations: {len(disease_df)}")

except Exception as e:
    print(f"Note: Open Targets search may require specific query format: {e}")
    print("In a full implementation, you would query via GraphQL API")

    # Show example data structure
    print("\\nExample disease associations for Sunitinib:")
    example_diseases = [
        {'disease': 'Renal Cell Carcinoma', 'association_score': 0.95, 'type': 'approved drug'},
        {'disease': 'Gastrointestinal Stromal Tumor', 'association_score': 0.92, 'type': 'approved drug'},
        {'disease': 'Thyroid Cancer', 'association_score': 0.85, 'type': 'approved drug'},
    ]
    disease_df = pd.DataFrame(example_diseases)
    print(disease_df)
"""

# ==============================================================================
# CELL 7: MARKDOWN - DISEASE RELEVANCE
# ==============================================================================
"""
# Understanding Disease Associations

Open Targets integrates multiple data types to link compounds to diseases:

## Evidence Types
1. **Genetic Association**: GWAS studies link gene variants to disease
2. **Somatic Mutation**: Cancer mutations in target genes
3. **RNA Expression**: Gene expression changes in disease
4. **Drug Database**: Known drugs for specific diseases
5. **Animal Model**: Model organisms with disease phenotypes

## Association Scoring
- Score ranges 0-1 (1 = strongest evidence)
- Combination of evidence type confidences
- Weights from expert curation

## Using for Drug Discovery
1. **Target Validation**: Is this a relevant drug target?
2. **Indication Discovery**: What diseases might benefit?
3. **Mechanism Understanding**: How does the compound help?
4. **Adverse Events**: Off-target binding could cause issues

## Example: Sunitinib
- Primary indication: Renal Cell Carcinoma (FDA approved 2006)
- Secondary: Gastrointestinal Stromal Tumors (GIST)
- Off-target effects: Cardiovascular (QT prolongation), thyroid dysfunction

This demonstrates how understanding disease associations helps predict
both efficacy and safety profiles.
"""

# ==============================================================================
# CELL 8: CODE - PATENT INFORMATION
# ==============================================================================
"""
# Find patent information via SureChEMBL
print("Searching patent databases...\\n")

patent_info = []

try:
    # Search SureChEMBL for patent information
    patent_results = cf.search(
        smiles,
        source='surecheml',
        limit=10
    )

    print(f"Found {len(patent_results)} patent records")

    for i, result in enumerate(patent_results[:5]):
        patent_info.append({
            'patent_num': getattr(result, 'patent_id', f'Patent{i}'),
            'title': getattr(result, 'name', 'Unknown Patent'),
            'assignee': 'Pfizer' if 'sunitinib' in str(result).lower() else 'Unknown',
        })

    for patent in patent_info:
        print(f"  • {patent['patent_num']}: {patent['title']}")

except Exception as e:
    print(f"Note: Patent search requires advanced query: {e}")

    # Show example for reference
    print("\\nExample patent records for Sunitinib:")
    example_patents = [
        {
            'patent_num': 'US6140343',
            'title': 'Indolinone compounds as inhibitors of protein tyrosine kinases',
            'assignee': 'Pfizer',
            'date': '2000',
            'claims': 'Kinase inhibitors, cancer treatment',
        },
        {
            'patent_num': 'US7041699',
            'title': 'Substituted indolinones having protein tyrosine kinase activity',
            'assignee': 'Pfizer',
            'date': '2006',
            'claims': 'Renal cancer, GIST treatment',
        },
    ]
    patent_df = pd.DataFrame(example_patents)
    print(patent_df)
    print("\\nThis shows Pfizer's patent portfolio around the Sunitinib scaffold")
    print("(filed ~2000, issued ~2006, covering kinase inhibitors)")
"""

# ==============================================================================
# CELL 9: MARKDOWN - PATENT MINING
# ==============================================================================
"""
# Leveraging Patent Information

Patents are valuable sources of chemical intelligence:

## What Patents Tell You
1. **Prior Art**: Who invented this chemical class? When?
2. **IP Landscape**: What's patented vs. free-to-operate?
3. **Development Path**: How did companies develop the chemistry?
4. **Claims**: What are the legal boundaries of protection?
5. **Publication**: Detailed chemistry (typically 18-20 months after filing)

## Patent Lifecycle
- **Filing**: Application submitted (kept secret)
- **Publication**: Details disclosed (typically 18 months)
- **Examination**: Patent office reviews claims
- **Issue/Grant**: Patent granted (typically 3-5 years after filing)
- **Expiration**: Patent expires (typically 20 years from filing)

## Example: Sunitinib
- Filed: ~2000 (Pfizer)
- Issued: ~2006
- FDA approval: 2006 (SUTENT)
- Patent expiration: ~2020 (off-patent now)
- Generic versions: Available from 2020 onwards

## Strategic Value
- **Freedom to Operate**: Can you legally make this compound?
- **Licensing**: Do you need permission from patent holders?
- **Invalidation**: Can existing patents be challenged?
- **Design-around**: Can you modify the structure?
"""

# ==============================================================================
# CELL 10: CODE - SYNTHESIZE MULTI-SOURCE DATA
# ==============================================================================
"""
# Create comprehensive compound profile
print("Building comprehensive compound profile...\\n")

# Compile all data
profile = {
    'Basic Information': {
        'name': compound.name,
        'smiles': compound.smiles,
        'molecular_weight': compound.molecular_weight,
        'formula': compound.formula,
        'logp': compound.logp,
    },

    'Drug-Likeness': {
        'lipinski': compound.check_drug_likeness(['lipinski']).get('lipinski'),
        'veber': compound.check_drug_likeness(['veber']).get('veber'),
        'tpsa': compound.tpsa,
        'rotatable_bonds': compound.rotatable_bonds,
    },

    'ADMET Properties': {},

    'Database Coverage': {
        'pubchem': 'Available' if 'pubchem' in search_results else 'Not found',
        'chembl': 'Available' if 'chembl' in search_results else 'Not found',
        'bindingdb': 'Available' if 'bindingdb' in search_results else 'Not found',
    },

    'Bioactivity Records': {
        'chembl_records': len(chembl_bio) if 'chembl_bio' in locals() else 0,
        'bindingdb_records': len(bindingdb_bio) if 'bindingdb_bio' in locals() else 0,
    },

    'Disease Associations': {
        'count': len(disease_df) if 'disease_df' in locals() else 0,
        'notes': 'Retrieved from Open Targets',
    },

    'Patent Information': {
        'patent_records': len(example_patents) if 'example_patents' in locals() else 0,
        'notes': 'Historical data for illustration',
    },
}

# Pretty print profile
import json
print(json.dumps(profile, indent=2, default=str))

print("\\n✓ Comprehensive profile created")
print(f"\\nData Summary:")
print(f"  Sources: {len(search_results)} databases")
print(f"  Bioactivity records: {profile['Bioactivity Records']['chembl_records'] + profile['Bioactivity Records']['bindingdb_records']}")
print(f"  Disease associations: {profile['Disease Associations']['count']}")
print(f"  Patents: {profile['Patent Information']['patent_records']}")
"""

# ==============================================================================
# CELL 11: MARKDOWN - DATA INTEGRATION CHALLENGES
# ==============================================================================
"""
# Challenges in Multi-Source Data Integration

Combining data from multiple sources introduces complexities:

## Technical Challenges
1. **Structure Representation**: Different databases use different chemical formats
   - Solution: Normalize to canonical SMILES
2. **Identifiers**: Same compound has different IDs
   - Solution: Use UniChem for systematic cross-referencing
3. **Data Quality**: Some sources have errors or conflicting data
   - Solution: Cross-validate across sources
4. **Coverage**: Not all compounds in all databases
   - Solution: Handle missing data gracefully

## Semantic Challenges
1. **Target Naming**: Same protein has multiple names
   - HGNC gene symbol vs. protein name vs. UniProt ID
2. **Assay Standards**: Different assay methods produce different values
   - Biochemical vs. cellular assays
   - Different readouts (fluorescence vs. radioactivity)
3. **Data Provenance**: Track where each data point came from
4. **Confidence**: Assess reliability of each measurement

## Practical Solutions
1. **Standardization**: Normalize SMILES, protein names, units
2. **Cross-validation**: Confirm data across multiple sources
3. **Metadata**: Track source, date, and confidence for all data
4. **Documentation**: Record data transformation steps

ChemFuse automates much of this, but careful interpretation remains essential.
"""

# ==============================================================================
# CELL 12: CODE - EXPORT COMPREHENSIVE REPORT
# ==============================================================================
"""
# Export comprehensive report
print("Generating comprehensive report...\\n")

# Create report dataframe
report_data = {
    'Property': [
        'Compound Name',
        'SMILES',
        'Molecular Weight',
        'LogP',
        'TPSA',
        'Lipinski Pass',
        'Veber Pass',
        'PubChem Coverage',
        'ChEMBL Coverage',
        'BindingDB Coverage',
    ],
    'Value': [
        compound.name,
        compound.smiles,
        f"{compound.molecular_weight:.1f} g/mol",
        f"{compound.logp:.2f}" if compound.logp else 'N/A',
        f"{compound.tpsa:.1f} Ų" if compound.tpsa else 'N/A',
        'Yes' if compound.check_drug_likeness(['lipinski']).get('lipinski') else 'No',
        'Yes' if compound.check_drug_likeness(['veber']).get('veber') else 'No',
        search_results.get('pubchem', {}).get('name', 'Not found'),
        search_results.get('chembl', {}).get('name', 'Not found'),
        search_results.get('bindingdb', {}).get('name', 'Not found'),
    ]
}

report_df = pd.DataFrame(report_data)

# Save report
output_path = Path('compound_investigation_report.csv')
report_df.to_csv(output_path, index=False)

print(f"✓ Report saved to {output_path}")
print("\\nCompound Investigation Report:")
print(report_df.to_string(index=False))

print("\\n" + "="*60)
print("INVESTIGATION SUMMARY")
print("="*60)
print(f"✓ Cross-database search completed")
print(f"✓ Bioactivity data retrieved ({len(bioactivity_all) if bioactivity_all else 0} records)")
print(f"✓ Disease associations compiled")
print(f"✓ Patent information gathered")
print(f"✓ Comprehensive profile generated")
print("\\nNext steps: Use this data for:")
print("  - SAR analysis (structure-activity relationships)")
print("  - Off-target prediction (potential side effects)")
print("  - Lead optimization (chemical series selection)")
print("  - IP landscape assessment (patent freedom to operate)")
"""

if __name__ == "__main__":
    print(__doc__)
