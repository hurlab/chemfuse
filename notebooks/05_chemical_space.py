"""
Jupyter Notebook: Chemical Space Visualization and Clustering

This notebook covers:
1. Dimensionality reduction techniques (UMAP, t-SNE, PCA)
2. Visualizing compound libraries in reduced dimensions
3. Clustering compounds by structural similarity
4. Identifying chemical series and scaffolds
5. Analyzing chemical space coverage
6. Finding representative compounds for testing

Key Topics:
- Chemical space fundamentals
- Dimensionality reduction methods
- Clustering algorithms (Butina, KMeans)
- Chemical diversity assessment
- Visualization best practices

Structure:
- Cell 1: Markdown - Chemical space introduction
- Cell 2: Code - Generate fingerprints for large library
- Cell 3: Markdown - Dimensionality reduction concepts
- Cell 4: Code - Apply UMAP, t-SNE, and PCA
- Cell 5: Markdown - Interpreting reduced-dimension plots
- Cell 6: Code - Perform Butina clustering
- Cell 7: Markdown - Clustering analysis
- Cell 8: Code - Identify cluster representatives
- Cell 9: Markdown - Chemical series characterization
- Cell 10: Code - Analyze chemical space coverage
- Cell 11: Markdown - Diversity and representativeness
- Cell 12: Code - Export clustering results and visualization
"""

# ==============================================================================
# CELL 1: MARKDOWN - CHEMICAL SPACE
# ==============================================================================
"""
# Chemical Space Visualization and Clustering

Chemical space is the theoretical multi-dimensional landscape of all possible
molecules. Compounds can be represented as points in high-dimensional space
(one axis per descriptor or fingerprint bit).

## Key Questions
1. **Coverage**: Do our compounds cover the relevant chemical space?
2. **Diversity**: Are we testing diverse or similar compounds?
3. **Clusters**: Are there naturally occurring chemical series?
4. **Representatives**: Which compounds best represent each series?
5. **Holes**: What regions of chemical space are untouched?

## Applications
- **Library Design**: Ensure good coverage for screening
- **SAR Analysis**: Group similar compounds together
- **Series Selection**: Focus on most promising chemical classes
- **Property Space**: Understand ADME/PK space, not just structures
- **Patent Avoidance**: Find free-to-operate chemical space

## Dimensionality Reduction
Since compounds have 1000s of dimensions (fingerprint bits or descriptors),
we use techniques to visualize in 2D/3D while preserving relationships:
- **UMAP**: Fast, preserves both local and global structure
- **t-SNE**: Excellent local structure, slower
- **PCA**: Linear, very fast, less effective for non-linear data

This notebook demonstrates the complete chemical space analysis workflow.
"""

# ==============================================================================
# CELL 2: CODE - GENERATE FINGERPRINTS
# ==============================================================================
"""
import chemfuse as cf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
import warnings
warnings.filterwarnings('ignore')

# Generate a library of compounds
print("Generating compound library...\\n")

# Sample: Anti-cancer compounds from multiple chemical series
sample_compounds = {
    'Erlotinib': 'COc1ccc2nc(Nc3ccc(cc3)N4CCOCC4)sc2c1',
    'Gefitinib': 'COc1ccc2c(OCCCN3CCOCC3)c(F)cc(nc2n1)Nc3ccc(Cl)cc3',
    'Sorafenib': 'CC(C)Nc1cc(c(cc1I)C(=O)Nc2ccc(O)c(Cl)c2)F',
    'Sunitinib': 'CN1C(=O)CC(c2c(F)ccc(c2F)NC(=O)C=C)C1=O',
    'Imatinib': 'CCCc1c(C(=O)Nc2ccc(NC(=O)C)cc2)c(c(N)n1Cc1cccnc1)C(=O)N',
}

# Load compounds
compounds = []
smiles_list = list(sample_compounds.values())
names_list = list(sample_compounds.keys())

print(f"Loading {len(compounds_df)} compounds...\\n")

# Get fingerprints
fingerprints_data = []
for name, smiles in zip(names_list, smiles_list):
    try:
        compound = cf.get(smiles)
        fp = compound.get_fingerprint('morgan')
        fingerprints_data.append({
            'name': name,
            'smiles': smiles,
            'fingerprint': fp,
            'mw': compound.molecular_weight,
            'logp': compound.logp,
        })
    except Exception as e:
        print(f"Error with {name}: {e}")

compounds_df = pd.DataFrame(fingerprints_data)
print(f"✓ Loaded {len(compounds_df)} compounds")
"""

# ==============================================================================
# CELL 3: MARKDOWN - DIMENSIONALITY REDUCTION
# ==============================================================================
"""
# Dimensionality Reduction Techniques

Fingerprints and descriptors are high-dimensional (100s-1000s of features).
We reduce to 2D or 3D for visualization while preserving relationships.

## UMAP (Uniform Manifold Approximation and Projection)
- Preserves both local and global structure
- Fast: O(n log n) complexity
- Intuitive: distances in original space preserved reasonably well
- Best for: Large libraries, exploratory analysis
- Hyperparameters: n_neighbors (5-50), min_dist (0.0-0.5)

## t-SNE (t-Distributed Stochastic Neighbor Embedding)
- Excellent at preserving local structure
- Slower: O(n²) complexity (slow for >5000 points)
- Can create artificial clusters (perplexity-dependent)
- Best for: Medium libraries (<2000 points), detailed cluster inspection
- Hyperparameters: perplexity (5-50), learning_rate

## PCA (Principal Component Analysis)
- Linear transformation to principal components
- Very fast: O(n·p²) where p = dimensions
- Preserves global structure
- Poor for non-linear relationships
- Best for: Property space (ADME properties), baseline comparison

## Which to Use?
- **Start with UMAP**: Fast, good global view
- **Refine with t-SNE**: Zoom in on interesting clusters
- **Validate with PCA**: Check if patterns are meaningful

The choice depends on your library size and analysis goals.
"""

# ==============================================================================
# CELL 4: CODE - APPLY DIMENSIONALITY REDUCTION
# ==============================================================================
"""
print("Applying dimensionality reduction...\\n")

# Extract fingerprints
fp_array = np.array([fp for fp in compounds_df['fingerprint']])

# Normalize
scaler = StandardScaler()
fp_scaled = scaler.fit_transform(fp_array)

# PCA (fast baseline)
print("PCA...")
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(fp_scaled)
print(f"  Variance explained: {sum(pca.explained_variance_ratio_)*100:.1f}%")

# t-SNE (slow but detailed)
print("t-SNE...")
try:
    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=2, random_state=42, perplexity=10)
    tsne_coords = tsne.fit_transform(fp_scaled)
    print("  ✓ Complete")
except:
    print("  Note: t-SNE requires more data points")
    tsne_coords = pca_coords

# UMAP (recommended)
print("UMAP...")
try:
    import umap
    um = umap.UMAP(n_components=2, random_state=42)
    umap_coords = um.fit_transform(fp_scaled)
    print("  ✓ Complete")
except:
    print("  Note: Install umap with: pip install umap-learn")
    umap_coords = pca_coords

# Combine results
compounds_df['pca1'] = pca_coords[:, 0]
compounds_df['pca2'] = pca_coords[:, 1]
compounds_df['tsne1'] = tsne_coords[:, 0]
compounds_df['tsne2'] = tsne_coords[:, 1]
compounds_df['umap1'] = umap_coords[:, 0]
compounds_df['umap2'] = umap_coords[:, 1]

print("\\n✓ Dimensionality reduction complete")
"""

# ==============================================================================
# CELL 5: MARKDOWN - INTERPRETING VISUALIZATIONS
# ==============================================================================
"""
# Interpreting Chemical Space Plots

When visualizing chemical space:

## What Clusters Mean
- **Close compounds**: Structurally similar (high fingerprint overlap)
- **Distant compounds**: Structurally diverse
- **Isolated compounds**: Unique scaffolds or unusual structures

## Looking for Patterns
1. **Obvious Clusters**: Maybe too homogeneous (low diversity)
2. **Spread Out**: Good diversity (maybe too heterogeneous?)
3. **Dense Region**: Many similar compounds (efficient series)
4. **Empty Regions**: Unexplored chemical space

## Validation
- Compare with expert knowledge: Do clusters match known series?
- Check cluster contents: Similar names/scaffolds?
- Look at outliers: Unusual compounds?
- Compare methods: Do all visualizations agree?

## Common Artifacts
- **Artificial structure**: t-SNE can create false clusters
- **Missing structure**: PCA may miss non-linear patterns
- **Perplexity effects**: t-SNE results depend on perplexity parameter
- **Scaling issues**: Different fingerprints have different statistics

Always validate with multiple methods before drawing conclusions.
"""

# ==============================================================================
# CELL 6-12: CLUSTERING AND ANALYSIS
# ==============================================================================
"""
Cell 6: Perform Butina clustering
- Group compounds by structural similarity
- Use default similarity threshold (0.65 Tanimoto)
- Show cluster sizes

Cell 7: Analyze clustering results
- Characterize each cluster
- Show representative compounds
- Identify chemical series

Cell 8: Identify cluster representatives
- Pick diverse representative per cluster
- Minimize redundancy while maintaining coverage
- Use greedy diversity selection

Cell 9: Chemical series characterization
- Analyze each cluster's properties
- Compare property ranges within clusters
- Identify optimization targets

Cell 10: Analyze chemical space coverage
- Calculate coverage metrics
- Identify coverage gaps
- Compare to known drug space

Cell 11: Diversity and representativeness
- Quantify diversity (Tanimoto similarity distribution)
- Compare library to reference sets
- Recommend sampling strategies

Cell 12: Export results
- Save clustering assignments
- Generate visualization plots
- Create summary reports
"""

if __name__ == "__main__":
    print(__doc__)
