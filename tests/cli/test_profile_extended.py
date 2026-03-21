"""Tests for profile command extended flags: --patents, --targets, --binding, --all."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest
from click.testing import CliRunner

from chemfuse.cli.main import cli
from chemfuse.models.compound import Compound, CompoundProperties
from chemfuse.models.patent import Patent
from chemfuse.models.target import TargetAssociation

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def aspirin_compound() -> Compound:
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        chembl_id="CHEMBL25",
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem"],
        properties=CompoundProperties(
            molecular_weight=180.16,
            xlogp=1.2,
            tpsa=63.6,
            hbd_count=1,
            hba_count=4,
        ),
    )


@pytest.fixture
def aspirin_with_patents() -> Compound:
    """Aspirin compound pre-loaded with patent data."""
    compound = Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        chembl_id="CHEMBL25",
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem", "surechembl"],
        patents=[
            Patent(
                patent_id="US-4568679-A",
                title="Coated aspirin tablet",
                filing_date="1984-03-15",
                assignee="Bayer AG",
                jurisdiction="US",
            )
        ],
    )
    return compound


@pytest.fixture
def aspirin_with_targets() -> Compound:
    """Aspirin compound pre-loaded with target association data."""
    compound = Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        chembl_id="CHEMBL25",
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem", "opentargets"],
        target_associations=[
            TargetAssociation(
                target_id="ENSG00000095303",
                target_name="PTGS1",
                disease_id="EFO_0003785",
                disease_name="pain",
                association_score=0.85,
                evidence_count=42,
            )
        ],
    )
    return compound


@pytest.fixture
def aspirin_enriched() -> Compound:
    """Aspirin compound with both patents and targets populated."""
    compound = Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        chembl_id="CHEMBL25",
        name="aspirin",
        formula="C9H8O4",
        sources=["pubchem", "surechembl", "opentargets"],
        patents=[
            Patent(
                patent_id="US-4568679-A",
                title="Coated aspirin tablet",
                filing_date="1984-03-15",
                assignee="Bayer AG",
                jurisdiction="US",
            )
        ],
        target_associations=[
            TargetAssociation(
                target_id="ENSG00000095303",
                target_name="PTGS1",
                disease_id="EFO_0003785",
                disease_name="pain",
                association_score=0.85,
                evidence_count=42,
            )
        ],
    )
    return compound


# ---------------------------------------------------------------------------
# Helper to capture enrich() kwargs
# ---------------------------------------------------------------------------


def _run_profile_with_enrich_capture(
    runner: CliRunner,
    compound: Compound,
    cli_args: list[str],
) -> tuple[object, dict]:
    """Invoke profile CLI and capture kwargs passed to compound.enrich().

    Returns (result, enrich_kwargs) where enrich_kwargs is the keyword
    arguments captured from the enrich() call, or {} if not called.
    """
    captured: dict = {}

    async def fake_enrich(self, *args, **kwargs):
        captured.update(kwargs)

    with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
        mock_get.return_value = compound
        with patch.object(Compound, "enrich", new=fake_enrich):
            result = runner.invoke(cli, cli_args)

    return result, captured


# ---------------------------------------------------------------------------
# Tests for --patents flag
# ---------------------------------------------------------------------------


class TestProfilePatentsFlag:
    """Tests for --patents flag on profile command."""

    def test_patents_flag_calls_enrich_with_patents(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--patents passes patents=True to compound.enrich()."""
        result, captured = _run_profile_with_enrich_capture(
            runner, aspirin_compound, ["profile", "2244", "--patents"]
        )
        assert result.exit_code == 0
        assert captured.get("patents") is True

    def test_patents_flag_exit_code_zero(
        self, runner: CliRunner, aspirin_with_patents: Compound
    ) -> None:
        """--patents flag produces exit code 0 with enriched compound."""
        result, _ = _run_profile_with_enrich_capture(
            runner, aspirin_with_patents, ["profile", "2244", "--patents"]
        )
        assert result.exit_code == 0

    def test_patents_flag_in_help(self, runner: CliRunner) -> None:
        """--patents flag appears in profile --help output."""
        result = runner.invoke(cli, ["profile", "--help"])
        assert result.exit_code == 0
        assert "--patents" in result.output

    def test_patents_section_in_plain_output(self) -> None:
        """_print_profile_plain includes Patents section when patents populated."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            name="aspirin",
            patents=[
                Patent(
                    patent_id="US-4568679-A",
                    title="Coated aspirin tablet",
                    filing_date="1984-03-15",
                    assignee="Bayer AG",
                )
            ],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "Patents" in combined
        assert "US-4568679-A" in combined

    def test_no_patents_section_when_empty(self) -> None:
        """_print_profile_plain omits Patents section when patents list is empty."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(smiles="CC", name="methane")

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "Patents" not in combined


# ---------------------------------------------------------------------------
# Tests for --targets flag
# ---------------------------------------------------------------------------


class TestProfileTargetsFlag:
    """Tests for --targets flag on profile command."""

    def test_targets_flag_calls_enrich_with_targets(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--targets passes targets=True to compound.enrich()."""
        result, captured = _run_profile_with_enrich_capture(
            runner, aspirin_compound, ["profile", "2244", "--targets"]
        )
        assert result.exit_code == 0
        assert captured.get("targets") is True

    def test_targets_flag_exit_code_zero(
        self, runner: CliRunner, aspirin_with_targets: Compound
    ) -> None:
        """--targets flag produces exit code 0."""
        result, _ = _run_profile_with_enrich_capture(
            runner, aspirin_with_targets, ["profile", "2244", "--targets"]
        )
        assert result.exit_code == 0

    def test_targets_flag_in_help(self, runner: CliRunner) -> None:
        """--targets flag appears in profile --help output."""
        result = runner.invoke(cli, ["profile", "--help"])
        assert result.exit_code == 0
        assert "--targets" in result.output

    def test_disease_associations_section_in_plain_output(self) -> None:
        """_print_profile_plain includes Disease Associations section when populated."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            name="aspirin",
            target_associations=[
                TargetAssociation(
                    target_id="ENSG00000095303",
                    target_name="PTGS1",
                    disease_name="pain",
                    association_score=0.85,
                )
            ],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "Disease Associations" in combined
        assert "PTGS1" in combined
        assert "pain" in combined

    def test_no_disease_section_when_empty(self) -> None:
        """_print_profile_plain omits Disease Associations section when list is empty."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(smiles="CC", name="methane")

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "Disease Associations" not in combined

    def test_association_score_in_plain_output(self) -> None:
        """_print_profile_plain shows association score in plain output."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            smiles="CC",
            name="methane",
            target_associations=[
                TargetAssociation(
                    target_id="T1",
                    target_name="GENE1",
                    disease_name="test disease",
                    association_score=0.75,
                )
            ],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "0.750" in combined


# ---------------------------------------------------------------------------
# Tests for --all flag
# ---------------------------------------------------------------------------


class TestProfileAllFlag:
    """Tests for --all flag enabling all enrichment sources."""

    def test_all_flag_calls_enrich_with_all_true(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--all passes patents=True, targets=True, binding=True to enrich()."""
        result, captured = _run_profile_with_enrich_capture(
            runner, aspirin_compound, ["profile", "2244", "--all"]
        )
        assert result.exit_code == 0
        assert captured.get("patents") is True
        assert captured.get("targets") is True
        assert captured.get("binding") is True

    def test_all_flag_exit_code_zero(
        self, runner: CliRunner, aspirin_enriched: Compound
    ) -> None:
        """--all flag produces exit code 0."""
        result, _ = _run_profile_with_enrich_capture(
            runner, aspirin_enriched, ["profile", "2244", "--all"]
        )
        assert result.exit_code == 0

    def test_all_flag_shown_in_help(self, runner: CliRunner) -> None:
        """--all flag appears in profile --help output."""
        result = runner.invoke(cli, ["profile", "--help"])
        assert result.exit_code == 0
        assert "--all" in result.output

    def test_all_flag_equivalent_to_individual_flags(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--all produces same enrich() kwargs as --patents --targets --binding."""
        _, kwargs_all = _run_profile_with_enrich_capture(
            runner, aspirin_compound, ["profile", "2244", "--all"]
        )
        _, kwargs_individual = _run_profile_with_enrich_capture(
            runner,
            aspirin_compound,
            ["profile", "2244", "--patents", "--targets", "--binding"],
        )
        assert kwargs_all.get("patents") == kwargs_individual.get("patents")
        assert kwargs_all.get("targets") == kwargs_individual.get("targets")
        assert kwargs_all.get("binding") == kwargs_individual.get("binding")


# ---------------------------------------------------------------------------
# Tests for --binding flag
# ---------------------------------------------------------------------------


class TestProfileBindingFlag:
    """Tests for --binding flag on profile command."""

    def test_binding_flag_calls_enrich_with_binding(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--binding passes binding=True to compound.enrich()."""
        result, captured = _run_profile_with_enrich_capture(
            runner, aspirin_compound, ["profile", "2244", "--binding"]
        )
        assert result.exit_code == 0
        assert captured.get("binding") is True

    def test_binding_flag_in_help(self, runner: CliRunner) -> None:
        """--binding flag appears in profile --help output."""
        result = runner.invoke(cli, ["profile", "--help"])
        assert result.exit_code == 0
        assert "--binding" in result.output


# ---------------------------------------------------------------------------
# Tests for combined flags
# ---------------------------------------------------------------------------


class TestProfileCombinedFlags:
    """Tests for combining multiple enrichment flags."""

    def test_patents_and_targets_combined(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """--patents --targets passes both flags to enrich()."""
        result, captured = _run_profile_with_enrich_capture(
            runner,
            aspirin_compound,
            ["profile", "2244", "--patents", "--targets"],
        )
        assert result.exit_code == 0
        assert captured.get("patents") is True
        assert captured.get("targets") is True

    def test_no_enrichment_flags_does_not_call_enrich(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """Without enrichment flags, enrich() is not called."""
        called = []

        async def spy_enrich(self, *args, **kwargs):
            called.append(kwargs)

        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            with patch.object(Compound, "enrich", new=spy_enrich):
                result = runner.invoke(cli, ["profile", "2244"])

        assert result.exit_code == 0
        assert len(called) == 0

    def test_enrichment_error_does_not_crash_profile(
        self, runner: CliRunner, aspirin_compound: Compound
    ) -> None:
        """Enrichment exception is caught; profile still exits 0."""

        async def failing_enrich(self, *args, **kwargs):
            raise Exception("Network error")

        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            with patch.object(Compound, "enrich", new=failing_enrich):
                result = runner.invoke(cli, ["profile", "2244", "--patents"])

        # Enrichment errors are caught as warnings, not fatal
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# Tests for plain output rendering with extended data
# ---------------------------------------------------------------------------


class TestProfilePlainOutputExtended:
    """Tests for _print_profile_plain rendering with patents and targets."""

    def test_plain_output_with_both_patents_and_targets(self) -> None:
        """_print_profile_plain renders both sections when both populated."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            name="aspirin",
            sources=["pubchem", "surechembl", "opentargets"],
            patents=[
                Patent(
                    patent_id="US-4568679-A",
                    title="Coated aspirin tablet",
                    filing_date="1984-03-15",
                )
            ],
            target_associations=[
                TargetAssociation(
                    target_id="ENSG00000095303",
                    target_name="PTGS1",
                    disease_name="pain",
                    association_score=0.85,
                )
            ],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "Patents" in combined
        assert "US-4568679-A" in combined
        assert "Disease Associations" in combined
        assert "PTGS1" in combined
        assert "pain" in combined

    def test_plain_output_none_association_score(self) -> None:
        """_print_profile_plain handles target association with None score."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            smiles="CC",
            name="methane",
            target_associations=[
                TargetAssociation(
                    target_id="T1",
                    target_name="GENE1",
                    disease_name="condition",
                    association_score=None,
                )
            ],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "N/A" in combined

    def test_plain_output_patent_missing_title(self) -> None:
        """_print_profile_plain handles patent with None title."""
        from chemfuse.cli.commands.profile import _print_profile_plain

        compound = Compound(
            smiles="CC",
            name="methane",
            patents=[Patent(patent_id="US-999", title=None, filing_date=None)],
        )

        output_lines = []

        def capture_echo(text="", **kwargs):
            output_lines.append(str(text))

        with patch("click.echo", side_effect=capture_echo):
            _print_profile_plain(compound)

        combined = "\n".join(output_lines)
        assert "US-999" in combined
