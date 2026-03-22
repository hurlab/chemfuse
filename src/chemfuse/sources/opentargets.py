"""Open Targets Platform GraphQL source adapter."""

from __future__ import annotations

import logging
from typing import Any

from chemfuse.core.exceptions import SourceError
from chemfuse.models.compound import Compound
from chemfuse.models.target import TargetAssociation
from chemfuse.sources._base import SourceAdapter
from chemfuse.sources._http import AsyncHTTPClient

OPENTARGETS_BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OPENTARGETS_RATE_LIMIT = 3.0  # conservative; no hard limit documented

logger = logging.getLogger(__name__)

# GraphQL query to retrieve disease associations for a drug by ChEMBL ID
_DRUG_DISEASE_QUERY = """
query DrugDiseaseAssociations($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    name
    linkedDiseases {
      count
      rows {
        disease {
          id
          name
        }
        score
        evidenceCount
      }
    }
    linkedTargets {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
      }
    }
  }
}
"""

# GraphQL query for target tractability by Ensembl gene ID
_TARGET_TRACTABILITY_QUERY = """
query TargetTractability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
  }
}
"""


class OpenTargetsAdapter(SourceAdapter):
    """Open Targets Platform GraphQL adapter.

    Retrieves disease-target associations, target tractability information,
    and genetic evidence. No API key required (CC0 license, public access).

    API: https://api.platform.opentargets.org/api/v4/graphql (POST)
    """

    name = "opentargets"
    base_url = OPENTARGETS_BASE_URL
    rate_limit = OPENTARGETS_RATE_LIMIT

    def __init__(self, cache: Any = None, timeout: float = 30.0) -> None:
        """Initialize the Open Targets adapter.

        Args:
            cache: Optional Cache instance.
            timeout: HTTP request timeout in seconds.
        """
        self._http = AsyncHTTPClient(
            source_name=self.name,
            base_url=self.base_url,
            rate_limit=self.rate_limit,
            timeout=timeout,
            cache=cache,
            headers={"Content-Type": "application/json"},
        )

    async def search(
        self,
        query: str,
        query_type: str = "name",
    ) -> list[Compound]:
        """Search is not the primary use case for Open Targets.

        Open Targets is used for target-disease association enrichment.
        Returns empty list.
        """
        return []

    async def get_by_id(self, identifier: str) -> Compound | None:
        """Not implemented for Open Targets.

        Returns:
            None always.
        """
        return None

    async def get_properties(self, identifier: str) -> dict:
        """Not implemented for Open Targets.

        Returns:
            Empty dict.
        """
        return {}

    def is_available(self) -> bool:
        """Always returns True."""
        return True

    async def search_by_chembl_id(
        self,
        chembl_id: str,
    ) -> list[TargetAssociation]:
        """Get disease-target associations for a drug by its ChEMBL ID.

        Queries the Open Targets Platform GraphQL API for all disease
        associations linked to the given drug.

        Args:
            chembl_id: ChEMBL compound identifier (e.g., "CHEMBL25").

        Returns:
            List of TargetAssociation objects, empty list if none found.

        Raises:
            SourceError: If the GraphQL response contains errors.
        """
        payload = {
            "query": _DRUG_DISEASE_QUERY,
            "variables": {"chemblId": chembl_id},
        }
        data = await self._http.post(self.base_url, data=payload)
        return self._parse_drug_disease_response(data, chembl_id)

    async def get_target_info(
        self,
        target_id: str,
    ) -> dict[str, Any]:
        """Get target tractability information for a given Ensembl gene ID.

        Args:
            target_id: Ensembl gene identifier (e.g., "ENSG00000095303").

        Returns:
            Dictionary with target details and tractability information.
        """
        payload = {
            "query": _TARGET_TRACTABILITY_QUERY,
            "variables": {"ensemblId": target_id},
        }
        data = await self._http.post(self.base_url, data=payload)
        # Check for GraphQL errors
        if isinstance(data, dict) and data.get("errors"):
            errors = data["errors"]
            msg = errors[0].get("message", "GraphQL error") if errors else "GraphQL error"
            raise SourceError(
                f"Open Targets GraphQL error: {msg}",
                source=self.name,
            )
        if isinstance(data, dict):
            return data.get("data", {}).get("target", {}) or {}
        return {}

    # --- Private helpers ---

    def _parse_drug_disease_response(
        self,
        data: Any,
        chembl_id: str,
    ) -> list[TargetAssociation]:
        """Parse Open Targets GraphQL drug/disease response into TargetAssociation objects.

        Args:
            data: Raw GraphQL response dict.
            chembl_id: The ChEMBL ID that was queried (for error context).

        Returns:
            List of TargetAssociation objects.

        Raises:
            SourceError: If the response contains GraphQL errors.
        """
        if not isinstance(data, dict):
            return []

        # R-OT-08: raise SourceError on GraphQL errors
        if data.get("errors"):
            errors = data["errors"]
            msg = errors[0].get("message", "GraphQL error") if errors else "GraphQL error"
            raise SourceError(
                f"Open Targets GraphQL error: {msg}",
                source=self.name,
            )

        drug_data = data.get("data", {})
        if not drug_data:
            return []

        drug = drug_data.get("drug")
        if not drug:
            # R-OT-07: no associations found -> return empty list
            return []

        associations: list[TargetAssociation] = []

        # Parse linkedDiseases which contains scored disease associations
        linked_diseases = drug.get("linkedDiseases") or {}
        disease_rows = linked_diseases.get("rows") or []

        # Build a target map from linkedTargets for lookup
        linked_targets = drug.get("linkedTargets") or {}
        target_rows = linked_targets.get("rows") or []
        # Map target_id -> target details for merging with diseases
        target_map: dict[str, dict[str, str]] = {}
        for tr in target_rows:
            t = tr.get("target") or {}
            tid = t.get("id", "")
            if tid:
                target_map[tid] = {
                    "id": tid,
                    "name": t.get("approvedSymbol") or t.get("approvedName") or tid,
                }

        if disease_rows:
            # When we have disease rows, create associations per disease
            # Use the first available target or a generic one
            # Most disease associations at the drug level don't have per-disease targets
            for row in disease_rows:
                if not isinstance(row, dict):
                    continue
                disease = row.get("disease") or {}
                disease_id = disease.get("id") or ""
                disease_name = disease.get("name") or ""
                score = row.get("score")
                evidence_count = row.get("evidenceCount")

                # Normalize score to [0, 1]
                normalized_score: float | None = None
                if score is not None:
                    try:
                        normalized_score = float(score)
                        normalized_score = max(0.0, min(1.0, normalized_score))
                    except (TypeError, ValueError):
                        normalized_score = None

                # If we have targets, create one entry per target-disease pair
                # Otherwise, create with a placeholder target
                if target_map:
                    for _tid, tinfo in target_map.items():
                        assoc = TargetAssociation(
                            target_id=tinfo["id"],
                            target_name=tinfo["name"],
                            disease_id=disease_id or None,
                            disease_name=disease_name or None,
                            association_score=normalized_score,
                            evidence_count=int(evidence_count) if evidence_count is not None else None,
                            source=self.name,
                        )
                        associations.append(assoc)
                else:
                    # No linked targets available; use drug name or chembl_id as placeholder
                    drug_name = drug.get("name") or chembl_id
                    assoc = TargetAssociation(
                        target_id=chembl_id,
                        target_name=drug_name,
                        disease_id=disease_id or None,
                        disease_name=disease_name or None,
                        association_score=normalized_score,
                        evidence_count=int(evidence_count) if evidence_count is not None else None,
                        source=self.name,
                    )
                    associations.append(assoc)

        elif target_rows:
            # Only targets available, no disease scores
            for tr in target_rows:
                t = tr.get("target") or {}
                tid = t.get("id", "")
                tsymbol = t.get("approvedSymbol") or t.get("approvedName") or tid
                if not tid:
                    continue
                assoc = TargetAssociation(
                    target_id=tid,
                    target_name=tsymbol,
                    source=self.name,
                )
                associations.append(assoc)

        return associations


__all__ = ["OpenTargetsAdapter"]
