"""
Jupyter Notebook: Screening Campaign

This notebook covers:
1. Loading compound data from CSV
2. Batch searching in databases
3. Predicting ADMET properties
4. Filtering by drug-likeness criteria
5. Clustering results by structure
6. Exporting prioritized compounds

Structure:
- Cell 1: Markdown - Overview of screening campaigns
- Cell 2: Code - Load sample compound library
- Cell 3: Markdown - Understanding the dataset
- Cell 4: Code - Batch search in PubChem
- Cell 5: Markdown - Enriching with properties
- Cell 6: Code - Predict ADMET properties
- Cell 7: Markdown - ADMET interpretation
- Cell 8: Code - Apply drug-likeness filters
- Cell 9: Markdown - Results filtering
- Cell 10: Code - Cluster by chemical structure
- Cell 11: Markdown - Interpreting clusters
- Cell 12: Code - Export results
"""

# ==============================================================================
# CELL 1: MARKDOWN - SCREENING CAMPAIGNS
# ==============================================================================
"""
# Screening Campaign: From Candidates to Leads

In this notebook, we'll conduct a virtual screening campaign to identify
promising drug candidates from a library of compounds.

## Workflow
1. Load compound candidates (e.g., hit compounds from an assay)
2. Retrieve detailed properties from public databases
3. Predict drug-like properties and ADMET
4. Filter against established criteria
5. Cluster by structure to identify chemical series
6. Prioritize candidates for synthesis/testing

## Typical Scenario
You have 100-1000 compounds from:
- High-throughput screening (HTS)
- Literature mining
- Vendor catalogs
- Computational design

Goal: Identify 10-20 promising candidates for further investigation.

This notebook demonstrates the complete workflow using ChemFuse.
"""

# ==============================================================================
# CELL 2: CODE - LOAD SAMPLE COMPOUNDS
# ==============================================================================
"""
import chemfuse as cf
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Load sample compound library
# In real scenario, you'd load from your experimental data
sample_smiles = [
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("Naproxen", "COc1ccc2cc(ccc2c1)C(C)C(=O)O"),
    ("Indomethacin", "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccccc1Cl"),
    ("Tolmetin", "CC(=O)Oc1ccccc1c1ccc(CC(=O)O)cc1"),
]

# Create dataframe
compounds_df = pd.DataFrame(sample_smiles, columns=['name', 'smiles'])
print(f"Loaded {len(compounds_df)} compounds")
print(compounds_df.head())

# In practice, you'd load from file:
# compounds_df = pd.read_csv('candidates.csv')
"""

# ==============================================================================
# CELL 3: MARKDOWN - UNDERSTANDING THE DATASET
# ==============================================================================
"""
# Understanding Your Compound Library

The input should be a CSV file with at least:
- **name** or **id**: Compound identifier
- **smiles**: SMILES string representation

Optional columns:
- **source**: Where the compound came from (HTS, literature, etc.)
- **potency**: Activity value if known
- **molecular_weight**: If pre-calculated
- **vendor**: Source vendor

In this example, we're using known NSAIDs (non-steroidal anti-inflammatory drugs)
as sample candidates. This allows us to validate our screening against known
good drugs.

The workflow will:
1. Enrich with comprehensive molecular properties
2. Predict ADMET characteristics
3. Apply drug-likeness filters
4. Identify chemical series by clustering
"""

# ==============================================================================
# CELL 4: CODE - BATCH SEARCH AND ENRICH
# ==============================================================================
"""
# Batch search in PubChem for enrichment
print("Searching PubChem for detailed properties...")
batch_results = cf.batch_search(
    compounds_df['smiles'].tolist(),
    source='pubchem',
    limit=1
)

print(f"Successfully matched {len(batch_results)} compounds")

# Enrich the original dataframe with PubChem data
enriched_data = []
for original_name, original_smiles in zip(compounds_df['name'], compounds_df['smiles']):
    try:
        compound = cf.get(original_smiles, source='pubchem')
        enriched_data.append({
            'input_name': original_name,
            'input_smiles': original_smiles,
            'pubchem_name': compound.name,
            'pubchem_id': compound.pubchem_cid,
            'smiles': compound.smiles,
            'formula': compound.formula,
            'mw': compound.molecular_weight,
            'logp': compound.logp,
        })
    except Exception as e:
        print(f"Error processing {original_name}: {e}")

enriched_df = pd.DataFrame(enriched_data)
print("\\nEnriched data:")
print(enriched_df[['input_name', 'pubchem_name', 'mw', 'logp']])
"""

# ==============================================================================
# CELL 5: MARKDOWN - PROPERTY ENRICHMENT
# ==============================================================================
"""
# Property Enrichment

By querying PubChem, we've enriched our dataset with:

- **Canonical SMILES**: Standardized representation
- **Molecular Weight**: Accurate calculation from structure
- **LogP**: Lipophilicity prediction
- **Molecular Formula**: Elemental composition

This standardization is important because:
1. Input SMILES may be non-canonical or contain errors
2. PubChem provides curated, validated structures
3. Standardized SMILES enables cross-database matching

Next, we'll calculate additional properties needed for drug-likeness
assessment and ADMET prediction.
"""

# ==============================================================================
# CELL 6: CODE - PREDICT ADMET PROPERTIES
# ==============================================================================
"""
# Predict ADMET properties for all compounds
print("Predicting ADMET properties...")

admet_results = []
for idx, row in enriched_df.iterrows():
    try:
        compound = cf.get(row['smiles'])

        # Predict ADMET
        admet = compound.predict_admet()

        # Get drug-likeness
        drug_like = compound.check_drug_likeness()

        # Get descriptors
        descriptors = compound.get_descriptors([
            'NumHBD', 'NumHBA', 'NumRotatableBonds', 'TPSA', 'NumAromaticRings'
        ])

        admet_results.append({
            'compound_idx': idx,
            'name': row['input_name'],
            # ADMET scores
            'absorption_score': admet.get('absorption', {}).get('score'),
            'distribution_score': admet.get('distribution', {}).get('score'),
            'metabolism_score': admet.get('metabolism', {}).get('score'),
            'excretion_score': admet.get('excretion', {}).get('score'),
            'toxicity_score': admet.get('toxicity', {}).get('score'),
            # Drug-likeness
            'lipinski_pass': drug_like.get('lipinski', False),
            'veber_pass': drug_like.get('veber', False),
            'ghose_pass': drug_like.get('ghose', False),
            # Descriptors
            'hbd': descriptors.get('NumHBD'),
            'hba': descriptors.get('NumHBA'),
            'rot_bonds': descriptors.get('NumRotatableBonds'),
            'tpsa': descriptors.get('TPSA'),
            'aromatic_rings': descriptors.get('NumAromaticRings'),
        })
    except Exception as e:
        print(f"Error with {row['input_name']}: {e}")

admet_df = pd.DataFrame(admet_results)
print("\\nADMET Predictions:")
print(admet_df[['name', 'absorption_score', 'distribution_score',
                 'lipinski_pass', 'veber_pass']])
"""

# ==============================================================================
# CELL 7: MARKDOWN - ADMET INTERPRETATION
# ==============================================================================
"""
# Understanding ADMET Predictions

**ADMET** stands for:
- **Absorption**: How well the compound is absorbed orally
- **Distribution**: How it spreads throughout the body
- **Metabolism**: How the body metabolizes it
- **Excretion**: How it's eliminated
- **Toxicity**: Any toxic effects

## Scoring
- Scores typically range from 0-1 (higher is better for ADME, lower for toxicity)
- Confidence intervals indicate prediction reliability
- Rule-based fallback used when ML model unavailable

## Drug-Likeness Filters
- **Lipinski Pass**: Passes Lipinski's Rule of Five (most important)
- **Veber Pass**: Stricter filter, good for oral bioavailability
- **Ghose Pass**: Optimization for druglikeness

## Typical Good Values
- All ADME scores > 0.6
- Lipinski and Veber passes = true
- LogP 0-5 range
- MW < 500
- TPSA 20-130
"""

# ==============================================================================
# CELL 8: CODE - APPLY FILTERS
# ==============================================================================
"""
# Merge ADMET with enriched data
screened_df = enriched_df.merge(
    admet_df,
    left_index=True,
    right_on='compound_idx'
)

# Apply filters
print("Applying drug-likeness filters...")

# Filter 1: Lipinski's Rule
pass_lipinski = screened_df[screened_df['lipinski_pass'] == True]
print(f"Pass Lipinski: {len(pass_lipinski)} / {len(screened_df)}")

# Filter 2: Multiple filters
pass_strict = pass_lipinski[
    (pass_lipinski['veber_pass'] == True) &
    (pass_lipinski['absorption_score'] > 0.5) &
    (pass_lipinski['toxicity_score'] > 0.6)
]
print(f"Pass Strict (Lipinski + Veber + ADMET): {len(pass_strict)} / {len(screened_df)}")

# Filter 3: Chemical space (MW, LogP ranges)
pass_chemical_space = pass_strict[
    (pass_strict['mw'] >= 200) & (pass_strict['mw'] <= 500) &
    (pass_strict['logp'] >= -1) & (pass_strict['logp'] <= 5)
]
print(f"Pass Chemical Space: {len(pass_chemical_space)} / {len(screened_df)}")

# Show passing compounds
print("\\nPassing compounds:")
print(pass_chemical_space[[
    'input_name', 'mw', 'logp', 'tpsa',
    'lipinski_pass', 'veber_pass', 'absorption_score'
]])
"""

# ==============================================================================
# CELL 9: MARKDOWN - FILTERING STRATEGY
# ==============================================================================
"""
# Filtering Strategy

The filtering approach above uses a staged pipeline:

## Stage 1: Basic Drug-Likeness (Lipinski)
- Removes obvious non-drug-like compounds
- ~90% typically pass

## Stage 2: Enhanced Filters
- Veber and ADMET thresholds
- Ensures good oral bioavailability
- ~70-80% of Lipinski pass

## Stage 3: Chemical Space
- Ensures similarity to known drugs
- Removes extreme LogP, MW values
- ~50-60% remain

## Key Metrics
- **Lipinski Pass**: Most important, quick initial filter
- **TPSA**: 20-130 ideal for good absorption
- **LogP**: Affects bioavailability and toxicity
- **Rotation Bonds**: Flexibility, <10 generally preferred

This staged approach efficiently prioritizes compounds while maintaining diversity.
"""

# ==============================================================================
# CELL 10: CODE - CLUSTER BY STRUCTURE
# ==============================================================================
"""
# Cluster passing compounds by structure
print("Clustering passing compounds...")

# Get compounds for clustering
passing_compounds = [
    cf.get(smiles) for smiles in pass_chemical_space['smiles']
]

if len(passing_compounds) > 1:
    # Perform clustering
    clusters = []
    cluster_results = []

    try:
        # Create a collection for clustering
        from chemfuse.core import CompoundCollection
        collection = CompoundCollection(passing_compounds)
        clusters = collection.cluster(method='butina')

        print(f"\\nIdentified {len(clusters)} chemical series:")

        for cluster_idx, cluster in enumerate(clusters):
            cluster_names = [c.name for c in cluster]
            print(f"\\nCluster {cluster_idx + 1}: {len(cluster)} compounds")
            print(f"  Representatives: {', '.join(cluster_names[:3])}")

            # Get cluster diversity (Tanimoto similarity)
            cluster_results.append({
                'cluster': cluster_idx,
                'size': len(cluster),
                'compounds': ', '.join(cluster_names)
            })
    except Exception as e:
        print(f"Clustering note: {e}")
        print("Using simple series grouping instead...")

        # Fallback: group by similar scaffolds
        for idx, compound in enumerate(passing_compounds):
            cluster_results.append({
                'cluster': idx,
                'size': 1,
                'compounds': compound.name
            })

cluster_summary = pd.DataFrame(cluster_results)
print("\\nCluster Summary:")
print(cluster_summary)
"""

# ==============================================================================
# CELL 11: MARKDOWN - INTERPRETING CLUSTERS
# ==============================================================================
"""
# Understanding Chemical Clusters

Clustering groups compounds by structural similarity. This is useful for:

## Why Clustering?
- **Series Development**: Identify which structural classes are most promising
- **Diversity**: Ensure candidates represent different mechanisms
- **SAR**: Build structure-activity relationships within clusters
- **Prioritization**: Focus on promising series

## Butina Clustering
- Distance-based algorithm
- Groups compounds by Tanimoto similarity of fingerprints
- Typical threshold: 0.65 (65% similar)
- Results in chemical series with related structures

## Common Series in Drug Discovery
- **Indole derivatives**: Common scaffold
- **Aryl acids**: Like NSAIDs in this example
- **Quinoline cores**: Antimalarials, kinase inhibitors
- **Imidazoles**: Antifungals, antihistamines

In this example, all compounds are aryl propionic acids (similar NSAID scaffolds),
so they cluster together. In a real screening campaign, you'd expect multiple
distinct chemical series.
"""

# ==============================================================================
# CELL 12: CODE - EXPORT PRIORITIZED RESULTS
# ==============================================================================
"""
# Create prioritized candidate list
print("Creating final candidate ranking...")

# Score compounds based on multiple criteria
def calculate_lead_score(row):
    \"\"\"Composite score for lead prioritization\"\"\"
    score = 0

    # Drug-likeness (max +30 points)
    if row['lipinski_pass']:
        score += 20
    if row['veber_pass']:
        score += 10

    # ADMET properties (max +40 points)
    admet_avg = np.mean([
        row['absorption_score'] or 0.5,
        row['distribution_score'] or 0.5,
        row['metabolism_score'] or 0.5,
        row['excretion_score'] or 0.5,
    ])
    score += admet_avg * 40

    # Avoid toxicity (max +20 points)
    if row['toxicity_score'] is not None and row['toxicity_score'] > 0.6:
        score += 20

    # Chemical space (max +10 points)
    if 0 <= row['logp'] <= 5:
        score += 5
    if 200 <= row['mw'] <= 500:
        score += 5

    return score

# Calculate scores
pass_chemical_space['lead_score'] = pass_chemical_space.apply(
    calculate_lead_score, axis=1
)

# Rank by score
ranked = pass_chemical_space.sort_values('lead_score', ascending=False)

# Export results
output_path = Path('screening_results.csv')
ranked[[
    'input_name', 'mw', 'logp', 'tpsa',
    'lipinski_pass', 'veber_pass',
    'absorption_score', 'toxicity_score',
    'lead_score'
]].to_csv(output_path, index=False)

print(f"\\n✓ Results saved to {output_path}")
print(f"\\nTop candidates (by lead score):")
print(ranked[[
    'input_name', 'mw', 'logp', 'lead_score'
]].head(10).to_string(index=False))

print(f"\\nScreening Summary:")
print(f"  Input compounds: {len(screened_df)}")
print(f"  Passed filters: {len(pass_chemical_space)}")
print(f"  Pass rate: {len(pass_chemical_space)/len(screened_df)*100:.1f}%")
"""

if __name__ == "__main__":
    print(__doc__)
