"""SDF file import and enriched export for ChemFuse."""
# @MX:ANCHOR: Public I/O boundary for SDF format; called by export.py and batch.py
# @MX:REASON: [AUTO] Fan-in >= 3: read_sdf and write_sdf are the canonical SDF API

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from chemfuse.models.compound import Compound, CompoundProperties

logger = logging.getLogger(__name__)

# Mapping of common SD tag names to CompoundProperties field names
_SD_PROP_MAP: dict[str, str] = {
    "MW": "molecular_weight",
    "MolecularWeight": "molecular_weight",
    "ExactMass": "exact_mass",
    "XLogP": "xlogp",
    "TPSA": "tpsa",
    "HBD": "hbd_count",
    "HBA": "hba_count",
    "RotatableBonds": "rotatable_bonds",
    "HeavyAtomCount": "heavy_atom_count",
    "AromaticRings": "aromatic_rings",
    "Complexity": "complexity",
    "SAScore": "sa_score",
    "fsp3": "fsp3",
}

# Mapping of CompoundProperties field names to SD tag names for writing
_PROP_TO_SD: dict[str, str] = {v: k for k, v in _SD_PROP_MAP.items() if v not in
                                 {_SD_PROP_MAP[k2] for k2 in _SD_PROP_MAP if k2 != k
                                  and _SD_PROP_MAP[k2] == v}}

# Use canonical SD tag names for writing
_PROP_TO_SD_CANONICAL: dict[str, str] = {
    "molecular_weight": "MW",
    "exact_mass": "ExactMass",
    "xlogp": "XLogP",
    "tpsa": "TPSA",
    "hbd_count": "HBD",
    "hba_count": "HBA",
    "rotatable_bonds": "RotatableBonds",
    "heavy_atom_count": "HeavyAtomCount",
    "aromatic_rings": "AromaticRings",
    "complexity": "Complexity",
    "sa_score": "SAScore",
    "fsp3": "fsp3",
}


def read_sdf(path: str) -> list[Compound]:
    """Read compounds from an SDF file.

    Parses each molecule in the SDF file, extracts its SMILES and any SD tag
    properties, and returns a list of Compound objects. Molecules that cannot
    be parsed by RDKit are skipped with a warning.

    Args:
        path: Path to the SDF file.

    Returns:
        List of Compound objects with SMILES and any SD properties populated.

    Raises:
        FileNotFoundError: If the path does not exist.
        OptionalDependencyError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem
    except ImportError as exc:
        from chemfuse.core.exceptions import OptionalDependencyError
        raise OptionalDependencyError("rdkit", "rdkit") from exc

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"SDF file not found: {path}")

    supplier = Chem.SDMolSupplier(str(p))
    compounds: list[Compound] = []
    skipped = 0

    for mol in supplier:
        if mol is None:
            skipped += 1
            continue

        smiles = Chem.MolToSmiles(mol)
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if name is not None:
            name = name.strip() or None

        # Extract all non-private SD properties
        raw_props: dict[str, Any] = {
            k: v for k, v in mol.GetPropsAsDict().items()
            if not k.startswith("_")
        }

        # Build compound kwargs
        compound_kwargs: dict[str, Any] = {
            "smiles": smiles,
            "sources": ["sdf"],
        }
        if name:
            compound_kwargs["name"] = name

        # Map well-known SD tags to identifier fields
        if "CID" in raw_props:
            try:
                compound_kwargs["cid"] = int(raw_props["CID"])
            except (ValueError, TypeError):
                pass

        if "ChEMBL_ID" in raw_props:
            compound_kwargs["chembl_id"] = str(raw_props["ChEMBL_ID"])

        if "InChIKey" in raw_props:
            compound_kwargs["inchikey"] = str(raw_props["InChIKey"])

        if "InChI" in raw_props:
            compound_kwargs["inchi"] = str(raw_props["InChI"])

        if "sources" in raw_props:
            extra_sources = [s.strip() for s in str(raw_props["sources"]).split(",") if s.strip()]
            compound_kwargs["sources"] = list(dict.fromkeys(["sdf"] + extra_sources))

        # Map SD property tags to CompoundProperties fields
        prop_kwargs: dict[str, Any] = {}
        for sd_tag, field_name in _SD_PROP_MAP.items():
            if sd_tag in raw_props and field_name not in prop_kwargs:
                try:
                    prop_kwargs[field_name] = float(raw_props[sd_tag])
                except (ValueError, TypeError):
                    pass

        compound_kwargs["properties"] = CompoundProperties(**prop_kwargs)

        # Preserve remaining SD tags as descriptors
        known_sd_keys = (
            {"CID", "ChEMBL_ID", "InChIKey", "InChI", "sources"}
            | set(_SD_PROP_MAP.keys())
        )
        extra_descriptors: dict[str, float] = {}
        for k, v in raw_props.items():
            if k not in known_sd_keys:
                if k.startswith("desc_"):
                    field = k[5:]
                    try:
                        extra_descriptors[field] = float(v)
                    except (ValueError, TypeError):
                        pass

        compound = Compound(**compound_kwargs)
        if extra_descriptors:
            compound.descriptors.update(extra_descriptors)

        compounds.append(compound)

    if skipped:
        logger.warning("Skipped %d unparseable molecules in %s", skipped, path)

    return compounds


def write_sdf(
    compounds: list[Compound],
    path: str,
    include_properties: bool = True,
) -> str:
    """Write compounds to an SDF file with all properties as SD tags.

    Each compound's SMILES is converted to a RDKit Mol object. When
    include_properties is True, all identifiers, physicochemical properties,
    computed descriptors, and source information are written as SD tags.

    Args:
        compounds: List of Compound objects to write.
        path: Output file path.
        include_properties: Whether to include all computed properties as SD tags.
            When False, only the molecule name is written.

    Returns:
        Absolute path to the written file as a string.

    Raises:
        OptionalDependencyError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem
    except ImportError as exc:
        from chemfuse.core.exceptions import OptionalDependencyError
        raise OptionalDependencyError("rdkit", "rdkit") from exc

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    writer = Chem.SDWriter(str(output_path))
    written = 0
    skipped = 0

    for compound in compounds:
        if not compound.smiles:
            skipped += 1
            continue

        mol = Chem.MolFromSmiles(compound.smiles)
        if mol is None:
            skipped += 1
            logger.warning("Could not parse SMILES for compound %r; skipping.", compound.name)
            continue

        # Set molecule name
        if compound.name:
            mol.SetProp("_Name", compound.name)

        if include_properties:
            # --- Identifiers ---
            if compound.cid is not None:
                mol.SetProp("CID", str(compound.cid))
            if compound.chembl_id:
                mol.SetProp("ChEMBL_ID", compound.chembl_id)
            if compound.inchikey:
                mol.SetProp("InChIKey", compound.inchikey)
            if compound.inchi:
                mol.SetProp("InChI", compound.inchi)

            # --- Physicochemical properties ---
            props = compound.properties.model_dump(exclude_none=True)
            for field_name, value in props.items():
                sd_tag = _PROP_TO_SD_CANONICAL.get(field_name, field_name)
                mol.SetProp(sd_tag, str(value))

            # --- Computed descriptors ---
            for desc_name, desc_value in compound.descriptors.items():
                mol.SetProp(f"desc_{desc_name}", str(desc_value))

            # --- Druglikeness summary ---
            if compound.druglikeness is not None:
                dl = compound.druglikeness
                overall = getattr(dl, "overall", None)
                if overall is not None:
                    mol.SetProp("druglikeness_pass", str(overall))

            # --- Source tracking ---
            if compound.sources:
                mol.SetProp("sources", ",".join(compound.sources))

        writer.write(mol)
        written += 1

    writer.close()

    if skipped:
        logger.warning("Skipped %d compounds with invalid/missing SMILES.", skipped)

    logger.debug("Wrote %d compounds to %s", written, output_path)
    return str(output_path.resolve())
