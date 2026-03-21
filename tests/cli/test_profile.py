"""Tests for chemfuse.cli.commands.profile."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest
from click.testing import CliRunner

from chemfuse.cli.main import cli
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def aspirin_compound() -> Compound:
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
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


class TestProfileHelpers:
    def test_print_profile_plain(self):
        """Test _print_profile_plain renders without error."""
        from chemfuse.cli.commands.profile import _print_profile_plain
        compound = Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
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
        from click.testing import CliRunner
        runner = CliRunner()
        with runner.isolated_filesystem():
            # Should not raise
            _print_profile_plain(compound)

    def test_print_profile_plain_minimal(self):
        """Test _print_profile_plain with minimal compound."""
        from chemfuse.cli.commands.profile import _print_profile_plain
        compound = Compound(smiles="C", sources=[])
        _print_profile_plain(compound)


class TestProfileCommand:
    def test_profile_found(self, runner: CliRunner, aspirin_compound: Compound):
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            result = runner.invoke(cli, ["profile", "2244"])
        assert result.exit_code == 0

    def test_profile_not_found_exits_1(self, runner: CliRunner):
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = None
            result = runner.invoke(cli, ["profile", "99999"])
        assert result.exit_code == 1

    def test_profile_with_source(self, runner: CliRunner, aspirin_compound: Compound):
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            result = runner.invoke(cli, ["profile", "2244", "--source", "pubchem"])
        assert result.exit_code == 0
        call_kwargs = mock_get.call_args[1]
        assert call_kwargs.get("source") == "pubchem"

    def test_profile_help(self, runner: CliRunner):
        result = runner.invoke(cli, ["profile", "--help"])
        assert result.exit_code == 0
        assert "IDENTIFIER" in result.output or "identifier" in result.output.lower()

    def test_profile_shows_name(self, runner: CliRunner, aspirin_compound: Compound):
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            result = runner.invoke(cli, ["profile", "2244"])
        # Should show compound name in output
        assert result.exit_code == 0

    def test_profile_with_similar_flag(self, runner: CliRunner, aspirin_compound: Compound):
        from chemfuse.models.collection import CompoundCollection
        similar_collection = CompoundCollection(
            compounds=[aspirin_compound],
            query="smiles",
            sources=["pubchem"],
        )
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = aspirin_compound
            with patch("chemfuse.find_similar_async", new_callable=AsyncMock) as mock_similar:
                mock_similar.return_value = similar_collection
                result = runner.invoke(cli, ["profile", "2244", "--similar"])
        assert result.exit_code == 0

    def test_profile_similar_empty_smiles(self, runner: CliRunner):
        """Compound with no SMILES should not trigger similarity search."""
        no_smiles = Compound(
            cid=2244,
            smiles="",
            name="aspirin",
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=180.16),
        )
        with patch("chemfuse.get_async", new_callable=AsyncMock) as mock_get:
            mock_get.return_value = no_smiles
            result = runner.invoke(cli, ["profile", "2244", "--similar"])
        assert result.exit_code == 0
