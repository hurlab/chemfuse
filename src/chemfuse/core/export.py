"""Export engine for CompoundCollection data."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from chemfuse.models.collection import CompoundCollection


def export_csv(collection: CompoundCollection, output_path: str | Path) -> Path:
    """Export a CompoundCollection to CSV format.

    Args:
        collection: CompoundCollection to export.
        output_path: Path for the output CSV file.

    Returns:
        Path to the created CSV file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = collection.to_dataframe()
    df.to_csv(output_path, index=False)
    return output_path


def export_json(
    collection: CompoundCollection,
    output_path: str | Path,
    indent: int = 2,
) -> Path:
    """Export a CompoundCollection to JSON format.

    Args:
        collection: CompoundCollection to export.
        output_path: Path for the output JSON file.
        indent: JSON indentation level.

    Returns:
        Path to the created JSON file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    for compound in collection:
        d = compound.to_dict()
        # Remove None values
        records.append({k: v for k, v in d.items() if v is not None})

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=indent, ensure_ascii=False)

    return output_path


def export_excel(
    collection: CompoundCollection,
    output_path: str | Path,
    sheet_name: str = "ChemFuse Results",
) -> Path:
    """Export a CompoundCollection to Excel format.

    Requires openpyxl to be installed.

    Args:
        collection: CompoundCollection to export.
        output_path: Path for the output Excel file.
        sheet_name: Name of the Excel sheet.

    Returns:
        Path to the created Excel file.

    Raises:
        ImportError: If openpyxl is not installed.
    """
    try:
        import openpyxl  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Excel export requires openpyxl. Install it with: pip install chemfuse[excel]"
        ) from exc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = collection.to_dataframe()
    df.to_excel(output_path, sheet_name=sheet_name, index=False)
    return output_path


def export_sdf(
    collection: CompoundCollection,
    output_path: str | Path,
    include_properties: bool = True,
) -> Path:
    """Export a CompoundCollection to an enriched SDF file.

    Delegates to :func:`chemfuse.core.sdf.write_sdf` so that all compound
    identifiers, physicochemical properties, computed descriptors, and source
    information are written as SD tags.

    Requires rdkit to be installed.

    Args:
        collection: CompoundCollection to export.
        output_path: Path for the output SDF file.
        include_properties: Whether to include all computed properties as SD
            tags.  Defaults to True (enriched output).

    Returns:
        Path to the created SDF file.

    Raises:
        OptionalDependencyError: If rdkit is not installed.
    """
    from chemfuse.core.sdf import write_sdf

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    write_sdf(list(collection), str(output_path), include_properties=include_properties)
    return output_path


def records_to_dataframe(records: list[dict[str, Any]]) -> pd.DataFrame:
    """Convert a list of record dicts to a DataFrame.

    Args:
        records: List of dictionaries.

    Returns:
        Pandas DataFrame.
    """
    if not records:
        return pd.DataFrame()
    df = pd.DataFrame(records)
    df = df.dropna(axis=1, how="all")
    return df
