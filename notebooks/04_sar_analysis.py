"""
Jupyter Notebook: SAR Analysis - Structure-Activity Relationships

This notebook covers:
1. Retrieving bioactivity data for a compound series
2. Analyzing structure-activity relationships (SAR)
3. Detecting activity cliffs (large potency changes from small structure changes)
4. Predicting activity based on molecular properties
5. Visualizing SAR patterns
6. Generating insights for medicinal chemistry optimization

Key Topics:
- SAR fundamentals (how structure affects activity)
- Activity cliffs (discontinuous SAR)
- Threshold effects and binding modes
- Optimization strategies based on SAR
- Case study: Kinase inhibitor series optimization

Structure:
- Cell 1: Markdown - SAR fundamentals
- Cell 2: Code - Load bioactivity data for a series
- Cell 3: Markdown - SAR patterns recognition
- Cell 4: Code - Detect activity cliffs
- Cell 5: Markdown - Interpreting activity cliffs
- Cell 6: Code - Analyze structural features vs activity
- Cell 7: Markdown - Feature importance
- Cell 8: Code - Predict and optimize
- Cell 9: Markdown - Optimization strategies
- Cell 10: Code - Generate SAR visualization
- Cell 11: Markdown - SAR application in lead optimization
- Cell 12: Code - Export SAR insights
"""

# ==============================================================================
# CELL 1: MARKDOWN - SAR FUNDAMENTALS
# ==============================================================================
"""
# Structure-Activity Relationships (SAR)

SAR analysis reveals how molecular structure determines biological activity.
This is the foundation of rational drug design.

## Key Concepts

### Activity Cliff
A dramatic change in potency from small structural modification.
- Example: Changing one atom can increase potency 100-fold
- Indicates different binding mode or critical feature
- Valuable for identifying optimization strategies

### SAR Rules
Patterns in how structure changes affect activity:
1. **Molecular weight**: Increasing MW often decreases activity (permeability issue)
2. **LogP**: Optimal range (usually 2-4); too high = toxicity, too low = low affinity
3. **TPSA**: Controls cell membrane permeability
4. **H-bonding**: Critical for target binding specificity
5. **Scaffolds**: Core structure determines target selectivity

### Medicinal Chemistry Optimization
1. **Hit to Lead**: Improve potency and selectivity
2. **Lead Optimization**: Improve ADMET properties
3. **Candidate Selection**: Find balance between all properties

## Why SAR Matters
- Guide synthesis (focus on promising analogs)
- Predict activity of untested compounds
- Identify failure modes (why some analogs are inactive)
- Design selectivity (target-specific binding)
"""

# ==============================================================================
# CELL 2: CODE - LOAD BIOACTIVITY DATA
# ==============================================================================
"""
import chemfuse as cf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# Load bioactivity data for a compound series
# Example: p38 kinase inhibitors from ChEMBL
print("Loading p38 kinase inhibitor bioactivity data...\\n")

# In a real scenario, you'd retrieve from ChEMBL via API
# For this example, we'll use sample data

sample_compounds = {
    'Compound1': {'smiles': 'CC(C)CC(C)(C)c1ccc(O)cc1', 'ic50_nm': 5, 'mw': 350, 'logp': 4.2},
    'Compound2': {'smiles': 'CC(C)CC(C)(C)c1ccc(OC)cc1', 'ic50_nm': 8, 'mw': 365, 'logp': 4.5},
    'Compound3': {'smiles': 'CC(C)CC(C)(C)c1ccc(OCC)cc1', 'ic50_nm': 12, 'mw': 380, 'logp': 4.8},
    'Compound4': {'smiles': 'CC(C)CC(C)(C)c1ccc(OCCC)cc1', 'ic50_nm': 25, 'mw': 395, 'logp': 5.1},
    'Compound5': {'smiles': 'CC(C)CC(C)c1ccc(OCCC)cc1', 'mw': 365, 'ic50_nm': 80, 'logp': 4.2},  # Activity cliff
}

bio_df = pd.DataFrame(sample_compounds).T
bio_df['log_ic50'] = np.log10(bio_df['ic50_nm'])
bio_df['potency'] = -np.log10(bio_df['ic50_nm'] / 1000)  # pIC50

print("Bioactivity Data (p38 Kinase Inhibitors):")
print(bio_df[['ic50_nm', 'potency', 'mw', 'logp']])
print(f"\\nCompound count: {len(bio_df)}")
print(f"Potency range: {bio_df['potency'].min():.2f} - {bio_df['potency'].max():.2f} pIC50")
"""

# ==============================================================================
# CELL 3: MARKDOWN - SAR PATTERNS
# ==============================================================================
"""
# Recognizing SAR Patterns

Common patterns in structure-activity data:

## 1. Linear SAR
Activity increases linearly with a property (e.g., MW, LogP)
- Easier to model and predict
- Less surprising activity changes
- Typical in early optimization

## 2. Parabolic SAR
Optimal range for a property
- Activity peaks at intermediate value
- Too low or too high = reduced activity
- Common for LogP, MW, TPSA

## 3. Activity Cliffs
Dramatic activity changes from small structural changes
- Indicates critical binding features
- Suggests different binding modes
- High value for optimization

## 4. Bimodal Distribution
Two distinct groups of actives
- May indicate different targets or mechanisms
- Separate optimization paths

## 5. Sparse Data
Incomplete structure-activity information
- May miss key SAR features
- Increases prediction uncertainty
"""

# ==============================================================================
# CELL 4: CODE - DETECT ACTIVITY CLIFFS
# ==============================================================================
"""
# Detect activity cliffs using structure similarity and activity difference
print("Detecting activity cliffs...\\n")

from sklearn.metrics.pairwise import cosine_similarity

# Calculate fingerprints for all compounds
fingerprints = []
for smiles in bio_df.index:
    try:
        compound = cf.get(smiles)
        fp = compound.get_fingerprint('morgan')
        fingerprints.append(fp)
    except:
        fingerprints.append(np.zeros(2048))

fingerprints = np.array(fingerprints)

# Calculate pairwise similarities
similarity_matrix = cosine_similarity(fingerprints)

# Find activity cliffs: high similarity but low activity correlation
cliffs = []
for i in range(len(bio_df)):
    for j in range(i+1, len(bio_df)):
        sim = similarity_matrix[i, j]
        activity_diff = abs(bio_df.iloc[i]['potency'] - bio_df.iloc[j]['potency'])

        # Activity cliff: similarity > 0.6 but activity difference > 1.0 pIC50
        if sim > 0.6 and activity_diff > 1.0:
            cliffs.append({
                'compound1': bio_df.index[i],
                'compound2': bio_df.index[j],
                'similarity': sim,
                'potency_diff': activity_diff,
                'potency1': bio_df.iloc[i]['potency'],
                'potency2': bio_df.iloc[j]['potency'],
            })

if cliffs:
    cliff_df = pd.DataFrame(cliffs)
    print(f"Detected {len(cliff_df)} activity cliffs:\\n")
    print(cliff_df[['compound1', 'compound2', 'similarity', 'potency_diff']])
else:
    print("No major activity cliffs detected in this series")
"""

# ==============================================================================
# CELL 5: MARKDOWN - INTERPRETING CLIFFS
# ==============================================================================
"""
# Understanding Activity Cliffs

Activity cliffs represent the most interesting SAR patterns.

## What They Mean
- Similar structures with very different potencies
- Indicates critical structural feature(s)
- May represent different binding modes
- High value for lead optimization

## Why They Occur
1. **Binding Mode Change**: Small change flips molecule orientation
2. **Critical Feature**: Loss of critical H-bond or interaction
3. **Allosteric Effect**: Change affects protein conformation
4. **Solubility Cliff**: Affects cellular permeability/exposure

## Using Cliffs for Optimization
1. Identify structural differences between high and low potency pairs
2. Test hypotheses about critical features
3. Design new analogs targeting identified features
4. Validate in biochemical assays

## Example: Kinase Inhibitors
- Losing a key H-bond to hinge region causes major potency loss
- Adding a hydrophobic substituent can increase potency 100-fold
- Increasing MW beyond threshold decreases cell permeability

Identifying and understanding these cliffs accelerates lead optimization.
"""

# ==============================================================================
# CELL 6-12: PLACEHOLDER CELLS
# ==============================================================================
"""
The remaining cells would cover:

Cell 6: Analyze structural features vs activity
- Correlation of MW, LogP, TPSA with potency
- Feature importance analysis
- Property ranges of active compounds

Cell 7: Feature importance and significance
- Statistical significance testing
- Which properties matter most?
- Target-specific requirements

Cell 8: Predict activity of new analogs
- Build SAR model
- Predict potency for designed compounds
- Confidence intervals

Cell 9: Optimization strategies
- Directional changes (increase/decrease MW, LogP)
- Scaffold modifications
- Series diversification

Cell 10: SAR visualization
- Scatter plots of key properties vs activity
- Structure-activity relationship plots
- Chemical space visualization

Cell 11: Application in lead optimization
- Decision making framework
- Risk assessment
- Synthesis prioritization

Cell 12: Export SAR insights
- SAR summary report
- Optimization recommendations
- Next synthesis targets
"""

if __name__ == "__main__":
    print(__doc__)
