"""Tests for chemfuse.cli._utils module."""

from __future__ import annotations

from pathlib import Path

import click
import pytest

from chemfuse.cli._utils import _validate_output_path


class TestValidateOutputPath:
    """Tests for _validate_output_path."""

    def test_accepts_normal_path(self) -> None:
        result = _validate_output_path(Path("output.csv"))
        assert result.name == "output.csv"

    def test_accepts_subdirectory_path(self) -> None:
        result = _validate_output_path(Path("results/output.csv"))
        assert result.name == "output.csv"

    def test_rejects_path_traversal(self) -> None:
        with pytest.raises(click.BadParameter, match="Path traversal not allowed"):
            _validate_output_path(Path("../../etc/passwd"))

    def test_rejects_embedded_traversal(self) -> None:
        with pytest.raises(click.BadParameter, match="Path traversal not allowed"):
            _validate_output_path(Path("results/../../secret.csv"))

    def test_returns_resolved_path(self) -> None:
        result = _validate_output_path(Path("output.csv"))
        assert result.is_absolute()
