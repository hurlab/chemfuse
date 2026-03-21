"""ChemFuse computation engine.

Provides local computation of molecular descriptors, fingerprints,
drug-likeness filters, QED scoring, and PAINS screening.

Key functions:
    compute_descriptors: Compute 200+ RDKit molecular descriptors.
    compute_descriptors_batch: Batch descriptor computation.
    compute_fingerprint: Generate molecular fingerprints.
    tanimoto_similarity: Compute Tanimoto similarity.
    bulk_tanimoto: Compute similarity against a list of targets.
    lipinski_filter: Lipinski Rule of Five.
    veber_filter: Veber oral bioavailability rules.
    ghose_filter: Ghose drug-likeness filter.
    egan_filter: Egan absorption filter.
    muegge_filter: Muegge drug-likeness filter.
    check_drug_likeness: Combined evaluation with all five filters.
    compute_qed: QED (Quantitative Estimate of Drug-likeness) score.
    screen_pains: PAINS substructure screening.
"""

from chemfuse.compute.descriptors import (
    compute_descriptors,
    compute_descriptors_batch,
    is_valid_smiles,
    smiles_to_inchi,
    smiles_to_inchikey,
)
from chemfuse.compute.druglikeness import (
    check_drug_likeness,
    egan_filter,
    ghose_filter,
    lipinski_filter,
    muegge_filter,
    veber_filter,
)
from chemfuse.compute.fingerprints import (
    bulk_tanimoto,
    compute_fingerprint,
    tanimoto_similarity,
)

__all__ = [
    # Descriptors
    "compute_descriptors",
    "compute_descriptors_batch",
    "is_valid_smiles",
    "smiles_to_inchi",
    "smiles_to_inchikey",
    # Fingerprints
    "compute_fingerprint",
    "tanimoto_similarity",
    "bulk_tanimoto",
    # Drug-likeness
    "lipinski_filter",
    "veber_filter",
    "ghose_filter",
    "egan_filter",
    "muegge_filter",
    "check_drug_likeness",
]
