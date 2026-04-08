"""SAR report generation for ChemFuse.

Generates publication-ready Structure-Activity Relationship reports
from compound collections with bioactivity data.
"""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.models.collection import CompoundCollection

logger = logging.getLogger(__name__)

# @MX:ANCHOR: generate_sar_report is the primary public API for SAR report generation
# @MX:REASON: Called by CompoundCollection.generate_sar_report and directly by users


def generate_sar_report(
    collection: CompoundCollection,
    activity_type: str = "IC50",
    target: str | None = None,
    format: str = "markdown",
) -> str:
    """Generate a SAR report from a compound collection.

    Analyzes scaffold distribution, activity trends, structure-activity
    relationships, and key molecular features across the collection.

    Args:
        collection: CompoundCollection with compounds (ideally with bioactivity data).
        activity_type: Activity type to analyze (IC50, Ki, etc.).
        target: Optional target name to filter bioactivities.
        format: Output format ("markdown" or "text").

    Returns:
        Formatted SAR report as a string.
    """
    n = len(collection)
    sources = collection.sources or []

    # --- Property statistics ---
    props_data = _collect_property_stats(collection)

    # --- Scaffold analysis ---
    scaffold_info = _collect_scaffold_info(collection)

    # --- Activity summary ---
    activity_info = _collect_activity_info(collection, activity_type, target)

    # --- Drug-likeness profile ---
    dl_profile = _collect_druglikeness_profile(collection)

    # --- Key findings ---
    findings = _generate_findings(
        n, scaffold_info, props_data, activity_info, dl_profile
    )

    if format == "markdown":
        return _render_markdown(
            n=n,
            sources=sources,
            activity_type=activity_type,
            target=target,
            scaffold_info=scaffold_info,
            props_data=props_data,
            activity_info=activity_info,
            dl_profile=dl_profile,
            findings=findings,
        )
    else:
        return _render_text(
            n=n,
            sources=sources,
            activity_type=activity_type,
            target=target,
            scaffold_info=scaffold_info,
            props_data=props_data,
            activity_info=activity_info,
            dl_profile=dl_profile,
            findings=findings,
        )


def _collect_property_stats(collection: CompoundCollection) -> dict[str, Any]:
    """Collect MW, logP, TPSA, HBD/HBA ranges and means from the collection."""
    mw_vals: list[float] = []
    logp_vals: list[float] = []
    tpsa_vals: list[float] = []
    hbd_vals: list[int] = []
    hba_vals: list[int] = []

    for compound in collection:
        p = compound.properties
        if p.molecular_weight is not None:
            mw_vals.append(p.molecular_weight)
        if p.xlogp is not None:
            logp_vals.append(p.xlogp)
        if p.tpsa is not None:
            tpsa_vals.append(p.tpsa)
        if p.hbd_count is not None:
            hbd_vals.append(p.hbd_count)
        if p.hba_count is not None:
            hba_vals.append(p.hba_count)

    def _stats(vals: list) -> dict[str, Any]:
        if not vals:
            return {"min": None, "max": None, "mean": None}
        return {
            "min": round(min(vals), 2),
            "max": round(max(vals), 2),
            "mean": round(sum(vals) / len(vals), 2),
        }

    return {
        "mw": _stats(mw_vals),
        "logp": _stats(logp_vals),
        "tpsa": _stats(tpsa_vals),
        "hbd": _stats(hbd_vals),
        "hba": _stats(hba_vals),
    }


def _collect_scaffold_info(collection: CompoundCollection) -> dict[str, Any]:
    """Collect scaffold frequency and diversity score."""
    top_scaffolds: list[dict] = []
    diversity: float | None = None
    unique_scaffolds: int = 0

    try:
        df = collection.scaffold_frequency()
        if not df.empty:
            unique_scaffolds = len(df)
            top_scaffolds = df.head(5).to_dict(orient="records")
    except Exception as exc:  # noqa: BLE001
        logger.debug("SAR report: scaffold analysis failed: %s", exc)

    try:
        diversity = collection.diversity_score()
    except Exception as exc:  # noqa: BLE001
        logger.debug("SAR report: diversity score failed: %s", exc)

    return {
        "unique_scaffolds": unique_scaffolds,
        "top_scaffolds": top_scaffolds,
        "diversity_score": diversity,
    }


def _collect_activity_info(
    collection: CompoundCollection,
    activity_type: str,
    target: str | None,
) -> dict[str, Any]:
    """Extract bioactivity summary from the collection."""
    activities: list[tuple[str, float]] = []  # (compound_name, value)

    for compound in collection:
        if not hasattr(compound, "bioactivities") or not compound.bioactivities:
            continue
        for bio in compound.bioactivities:
            bio_type = getattr(bio, "activity_type", None) or getattr(bio, "type", None) or ""
            if activity_type.lower() not in str(bio_type).lower():
                continue
            # Optional target filter
            if target is not None:
                bio_target = getattr(bio, "target_name", None) or getattr(bio, "target", None) or ""
                if target.lower() not in str(bio_target).lower():
                    continue
            value = getattr(bio, "value", None) or getattr(bio, "value_nm", None)
            if value is not None:
                try:
                    name = compound.name or compound.smiles or "unknown"
                    activities.append((name, float(value)))
                except (TypeError, ValueError):
                    pass
            break  # Use first matching bioactivity per compound

    if not activities:
        return {"has_data": False}

    values = [v for _, v in activities]
    sorted_by_activity = sorted(activities, key=lambda x: x[1])

    return {
        "has_data": True,
        "count": len(activities),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
        "mean": round(sum(values) / len(values), 4),
        "most_potent": sorted_by_activity[:3],
        "least_potent": sorted_by_activity[-3:],
    }


def _collect_druglikeness_profile(collection: CompoundCollection) -> dict[str, Any]:
    """Collect percentage of compounds passing each drug-likeness filter."""
    if not collection.compounds:
        return {}

    counts = {"lipinski": 0, "veber": 0, "ghose": 0, "egan": 0, "muegge": 0}
    evaluated = 0

    for compound in collection:
        dl = compound.druglikeness
        if dl is None:
            try:
                compound.check_drug_likeness()
                dl = compound.druglikeness
            except Exception as exc:  # noqa: BLE001
                logger.debug("SAR report: drug-likeness failed for %s: %s", compound.name, exc)
                continue

        if dl is None:
            continue

        evaluated += 1
        if dl.lipinski.pass_filter:
            counts["lipinski"] += 1
        if dl.veber.pass_filter:
            counts["veber"] += 1
        if dl.ghose.pass_filter:
            counts["ghose"] += 1
        if dl.egan.pass_filter:
            counts["egan"] += 1
        if dl.muegge.pass_filter:
            counts["muegge"] += 1

    if evaluated == 0:
        return {"evaluated": 0}

    return {
        "evaluated": evaluated,
        "lipinski_pct": round(counts["lipinski"] / evaluated * 100, 1),
        "veber_pct": round(counts["veber"] / evaluated * 100, 1),
        "ghose_pct": round(counts["ghose"] / evaluated * 100, 1),
        "egan_pct": round(counts["egan"] / evaluated * 100, 1),
        "muegge_pct": round(counts["muegge"] / evaluated * 100, 1),
    }


def _generate_findings(
    n: int,
    scaffold_info: dict,
    props_data: dict,
    activity_info: dict,
    dl_profile: dict,
) -> list[str]:
    """Auto-generate bullet point findings from the report data."""
    findings: list[str] = []

    if n == 0:
        findings.append("The collection is empty.")
        return findings

    # Scaffold diversity
    unique = scaffold_info.get("unique_scaffolds", 0)
    if unique > 0 and n > 0:
        ratio = round(unique / n * 100, 1)
        if ratio >= 80:
            findings.append(
                f"High scaffold diversity: {unique} unique scaffolds among {n} compounds ({ratio}%)."
            )
        elif ratio <= 20:
            findings.append(
                f"Low scaffold diversity: only {unique} unique scaffolds among {n} compounds ({ratio}%), "
                "suggesting a focused, scaffold-based series."
            )
        else:
            findings.append(
                f"Moderate scaffold diversity: {unique} unique scaffolds across {n} compounds ({ratio}%)."
            )

    # Diversity score
    div_score = scaffold_info.get("diversity_score")
    if div_score is not None:
        if div_score >= 0.7:
            findings.append(f"High chemical diversity (mean Tanimoto distance: {div_score:.2f}).")
        elif div_score <= 0.3:
            findings.append(
                f"Low chemical diversity (mean Tanimoto distance: {div_score:.2f}), "
                "consistent with a congeneric series."
            )

    # Molecular weight
    mw = props_data.get("mw", {})
    if mw.get("mean") is not None:
        if mw["mean"] > 500:
            findings.append(
                f"Mean MW ({mw['mean']} Da) exceeds Lipinski's 500 Da threshold; "
                "oral bioavailability may be a concern."
            )
        else:
            findings.append(
                f"Mean MW ({mw['mean']} Da) is within Lipinski's 500 Da threshold."
            )

    # LogP
    logp = props_data.get("logp", {})
    if logp.get("max") is not None and logp["max"] > 5:
        findings.append(
            f"Maximum logP ({logp['max']}) exceeds 5; some compounds may have poor solubility."
        )

    # Activity
    if activity_info.get("has_data"):
        most_potent = activity_info["most_potent"]
        if most_potent:
            name, val = most_potent[0]
            findings.append(
                f"Most potent compound: {name} ({val} nM)."
            )
        findings.append(
            f"Activity range: {activity_info['min']} – {activity_info['max']} nM "
            f"(mean: {activity_info['mean']} nM)."
        )
    else:
        findings.append("No bioactivity data found in the collection for the specified activity type.")

    # Drug-likeness
    if dl_profile.get("evaluated", 0) > 0:
        lip_pct = dl_profile.get("lipinski_pct", 0)
        if lip_pct == 100:
            findings.append("All compounds pass Lipinski's Rule of Five.")
        elif lip_pct >= 80:
            findings.append(f"{lip_pct}% of compounds pass Lipinski's Rule of Five.")
        elif lip_pct < 50:
            findings.append(
                f"Only {lip_pct}% of compounds pass Lipinski's Rule of Five; "
                "drug-likeness may require optimization."
            )

    return findings


def _render_markdown(
    n: int,
    sources: list[str],
    activity_type: str,
    target: str | None,
    scaffold_info: dict,
    props_data: dict,
    activity_info: dict,
    dl_profile: dict,
    findings: list[str],
) -> str:
    """Render the SAR report in Markdown format."""
    lines: list[str] = []

    # Title
    target_str = f" — Target: {target}" if target else ""
    lines.append(f"# SAR Report{target_str}")
    lines.append("")

    # Executive Summary
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(f"- **Compound count**: {n}")
    lines.append(f"- **Source databases**: {', '.join(sources) if sources else 'N/A'}")
    lines.append(f"- **Activity type analyzed**: {activity_type}")
    unique = scaffold_info.get("unique_scaffolds", 0)
    lines.append(f"- **Unique scaffolds**: {unique}")
    if activity_info.get("has_data"):
        lines.append(
            f"- **Activity range**: {activity_info['min']} – {activity_info['max']} nM"
        )
    else:
        lines.append("- **Activity range**: N/A (no bioactivity data)")
    lines.append("")

    # Scaffold Analysis
    lines.append("## Scaffold Analysis")
    lines.append("")
    div = scaffold_info.get("diversity_score")
    div_str = f"{div:.3f}" if div is not None else "N/A"
    lines.append(f"- **Scaffold diversity score** (mean Tanimoto distance): {div_str}")
    lines.append(f"- **Unique scaffolds**: {unique}")
    lines.append("")
    top = scaffold_info.get("top_scaffolds", [])
    if top:
        lines.append("### Top Scaffolds by Frequency")
        lines.append("")
        lines.append("| Scaffold | Count | Percentage |")
        lines.append("|----------|-------|------------|")
        for row in top:
            lines.append(
                f"| `{row.get('scaffold', '')}` | {row.get('count', '')} | {row.get('percentage', '')}% |"
            )
        lines.append("")

    # Property Distribution
    lines.append("## Property Distribution")
    lines.append("")

    def _fmt_stats(stats: dict) -> str:
        if stats.get("mean") is None:
            return "N/A"
        return f"min {stats['min']}, mean {stats['mean']}, max {stats['max']}"

    lines.append(f"- **MW (Da)**: {_fmt_stats(props_data.get('mw', {}))}")
    lines.append(f"- **logP**: {_fmt_stats(props_data.get('logp', {}))}")
    lines.append(f"- **TPSA (Å²)**: {_fmt_stats(props_data.get('tpsa', {}))}")
    lines.append(f"- **HBD**: {_fmt_stats(props_data.get('hbd', {}))}")
    lines.append(f"- **HBA**: {_fmt_stats(props_data.get('hba', {}))}")
    lines.append("")

    # Activity Summary
    lines.append("## Activity Summary")
    lines.append("")
    if activity_info.get("has_data"):
        lines.append(f"- **Compounds with {activity_type} data**: {activity_info['count']}")
        lines.append(f"- **Range**: {activity_info['min']} – {activity_info['max']} nM")
        lines.append(f"- **Mean**: {activity_info['mean']} nM")
        lines.append("")
        most_potent = activity_info.get("most_potent", [])
        if most_potent:
            lines.append("### Most Potent Compounds")
            lines.append("")
            for name, val in most_potent:
                lines.append(f"- {name}: {val} nM")
            lines.append("")
        least_potent = activity_info.get("least_potent", [])
        if least_potent:
            lines.append("### Least Potent Compounds")
            lines.append("")
            for name, val in least_potent:
                lines.append(f"- {name}: {val} nM")
            lines.append("")
    else:
        lines.append(f"No {activity_type} bioactivity data found in this collection.")
        lines.append("")

    # Drug-likeness Profile
    lines.append("## Drug-likeness Profile")
    lines.append("")
    if dl_profile.get("evaluated", 0) > 0:
        lines.append(f"- **Evaluated**: {dl_profile['evaluated']} compounds")
        lines.append(f"- **Lipinski (Ro5)**: {dl_profile.get('lipinski_pct', 0)}% pass")
        lines.append(f"- **Veber**: {dl_profile.get('veber_pct', 0)}% pass")
        lines.append(f"- **Ghose**: {dl_profile.get('ghose_pct', 0)}% pass")
        lines.append(f"- **Egan**: {dl_profile.get('egan_pct', 0)}% pass")
        lines.append(f"- **Muegge**: {dl_profile.get('muegge_pct', 0)}% pass")
    else:
        lines.append("Drug-likeness data not available (RDKit required).")
    lines.append("")

    # Key Findings
    lines.append("## Key Findings")
    lines.append("")
    for finding in findings:
        lines.append(f"- {finding}")
    lines.append("")

    return "\n".join(lines)


def _render_text(
    n: int,
    sources: list[str],
    activity_type: str,
    target: str | None,
    scaffold_info: dict,
    props_data: dict,
    activity_info: dict,
    dl_profile: dict,
    findings: list[str],
) -> str:
    """Render the SAR report in plain text format."""
    lines: list[str] = []
    sep = "=" * 60

    target_str = f" - Target: {target}" if target else ""
    lines.append(f"SAR REPORT{target_str}")
    lines.append(sep)
    lines.append("")

    lines.append("EXECUTIVE SUMMARY")
    lines.append("-" * 40)
    lines.append(f"Compound count:       {n}")
    lines.append(f"Source databases:     {', '.join(sources) if sources else 'N/A'}")
    lines.append(f"Activity type:        {activity_type}")
    lines.append(f"Unique scaffolds:     {scaffold_info.get('unique_scaffolds', 0)}")
    if activity_info.get("has_data"):
        lines.append(
            f"Activity range:       {activity_info['min']} - {activity_info['max']} nM"
        )
    else:
        lines.append("Activity range:       N/A")
    lines.append("")

    lines.append("SCAFFOLD ANALYSIS")
    lines.append("-" * 40)
    div = scaffold_info.get("diversity_score")
    lines.append(f"Diversity score:      {f'{div:.3f}' if div is not None else 'N/A'}")
    lines.append(f"Unique scaffolds:     {scaffold_info.get('unique_scaffolds', 0)}")
    lines.append("")

    lines.append("PROPERTY DISTRIBUTION")
    lines.append("-" * 40)

    def _fmt(stats: dict) -> str:
        if stats.get("mean") is None:
            return "N/A"
        return f"min={stats['min']}, mean={stats['mean']}, max={stats['max']}"

    lines.append(f"MW (Da):    {_fmt(props_data.get('mw', {}))}")
    lines.append(f"logP:       {_fmt(props_data.get('logp', {}))}")
    lines.append(f"TPSA:       {_fmt(props_data.get('tpsa', {}))}")
    lines.append(f"HBD:        {_fmt(props_data.get('hbd', {}))}")
    lines.append(f"HBA:        {_fmt(props_data.get('hba', {}))}")
    lines.append("")

    lines.append("ACTIVITY SUMMARY")
    lines.append("-" * 40)
    if activity_info.get("has_data"):
        lines.append(f"Count: {activity_info['count']}")
        lines.append(f"Range: {activity_info['min']} - {activity_info['max']} nM")
        lines.append(f"Mean:  {activity_info['mean']} nM")
    else:
        lines.append(f"No {activity_type} data available.")
    lines.append("")

    lines.append("DRUG-LIKENESS PROFILE")
    lines.append("-" * 40)
    if dl_profile.get("evaluated", 0) > 0:
        lines.append(f"Lipinski: {dl_profile.get('lipinski_pct', 0)}% pass")
        lines.append(f"Veber:    {dl_profile.get('veber_pct', 0)}% pass")
        lines.append(f"Ghose:    {dl_profile.get('ghose_pct', 0)}% pass")
        lines.append(f"Egan:     {dl_profile.get('egan_pct', 0)}% pass")
        lines.append(f"Muegge:   {dl_profile.get('muegge_pct', 0)}% pass")
    else:
        lines.append("Not available.")
    lines.append("")

    lines.append("KEY FINDINGS")
    lines.append("-" * 40)
    for finding in findings:
        lines.append(f"* {finding}")
    lines.append("")

    return "\n".join(lines)
