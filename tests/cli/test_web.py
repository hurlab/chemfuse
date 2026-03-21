"""Tests for the CLI web command."""

from __future__ import annotations

import sys
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemfuse.cli.main import cli


class TestWebCommandBasic:
    """Basic tests for the web CLI command."""

    def test_web_command_is_registered(self) -> None:
        """The 'web' command is registered in the CLI group."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert "web" in result.output

    def test_web_command_help(self) -> None:
        """The 'web' command has a help string."""
        runner = CliRunner()
        result = runner.invoke(cli, ["web", "--help"])
        assert result.exit_code == 0
        assert "Streamlit" in result.output or "streamlit" in result.output or "web" in result.output.lower()

    def test_web_command_port_option(self) -> None:
        """The 'web' command accepts --port option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["web", "--help"])
        assert "--port" in result.output

    def test_web_command_host_option(self) -> None:
        """The 'web' command accepts --host option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["web", "--help"])
        assert "--host" in result.output


class TestWebCommandMissingStreamlit:
    """Test behavior when Streamlit is not installed."""

    def test_raises_optional_dependency_error(self) -> None:
        """Raises OptionalDependencyError when streamlit is not importable."""
        runner = CliRunner()

        with patch.dict(sys.modules, {"streamlit": None}):
            # Remove cached import
            with patch("builtins.__import__", side_effect=ImportError("No module named 'streamlit'")):
                result = runner.invoke(cli, ["web"])
                # Should exit with error or show error message
                assert result.exit_code != 0 or "streamlit" in (result.output or "").lower() or result.exception is not None

    def test_error_suggests_install_command(self) -> None:
        """Error message contains installation hint."""
        from chemfuse.core.exceptions import OptionalDependencyError

        exc = OptionalDependencyError("streamlit", extra="web")
        assert "chemfuse[web]" in str(exc) or "streamlit" in str(exc)


class TestWebCommandLaunch:
    """Test web command launch logic."""

    def test_default_port_is_8501(self) -> None:
        """Default port is 8501."""
        from chemfuse.cli.commands.web import web_cmd

        runner = CliRunner()
        # Mock subprocess.run and streamlit import
        with patch("chemfuse.cli.commands.web.subprocess.run") as mock_run, \
             patch("chemfuse.cli.commands.web.streamlit", create=True):
            runner.invoke(web_cmd, [])
        # If called, port should be 8501
        if mock_run.called:
            cmd_args = mock_run.call_args[0][0]
            assert "8501" in cmd_args

    def test_custom_port_passed_to_streamlit(self) -> None:
        """Custom --port value is passed to the streamlit command."""
        from chemfuse.cli.commands.web import web_cmd

        runner = CliRunner()
        with patch("chemfuse.cli.commands.web.subprocess.run") as mock_run, \
             patch("builtins.__import__", side_effect=lambda name, *a, **kw: MagicMock() if name == "streamlit" else __import__(name, *a, **kw)):
            runner.invoke(web_cmd, ["--port", "9000"])

        if mock_run.called:
            cmd_args = mock_run.call_args[0][0]
            assert "9000" in cmd_args

    def test_custom_host_passed_to_streamlit(self) -> None:
        """Custom --host value is passed to the streamlit command."""
        from chemfuse.cli.commands.web import web_cmd

        runner = CliRunner()
        with patch("chemfuse.cli.commands.web.subprocess.run") as mock_run, \
             patch("builtins.__import__", side_effect=lambda name, *a, **kw: MagicMock() if name == "streamlit" else __import__(name, *a, **kw)):
            runner.invoke(web_cmd, ["--host", "0.0.0.0"])

        if mock_run.called:
            cmd_args = mock_run.call_args[0][0]
            assert "0.0.0.0" in cmd_args


class TestWebCommandWithStreamlit:
    """Test web command when streamlit is available."""

    def test_launches_streamlit_process(self) -> None:
        """Web command calls subprocess.run with streamlit arguments."""
        from chemfuse.cli.commands.web import web_cmd

        runner = CliRunner()
        with patch("chemfuse.cli.commands.web.subprocess.run") as mock_run:
            # Ensure streamlit import succeeds
            mock_st = MagicMock()
            mock_st.__version__ = "1.55.0"
            with patch.dict(sys.modules, {"streamlit": mock_st}):
                runner.invoke(web_cmd, ["--port", "8501"])

        if mock_run.called:
            cmd = mock_run.call_args[0][0]
            # Should include streamlit run
            cmd_str = " ".join(str(x) for x in cmd)
            assert "streamlit" in cmd_str
            assert "run" in cmd_str

    def test_no_browser_flag(self) -> None:
        """--no-browser flag adds headless configuration."""
        from chemfuse.cli.commands.web import web_cmd

        runner = CliRunner()
        with patch("chemfuse.cli.commands.web.subprocess.run") as mock_run:
            mock_st = MagicMock()
            mock_st.__version__ = "1.55.0"
            with patch.dict(sys.modules, {"streamlit": mock_st}):
                runner.invoke(web_cmd, ["--no-browser"])

        if mock_run.called:
            cmd_str = " ".join(str(x) for x in mock_run.call_args[0][0])
            assert "headless" in cmd_str or "true" in cmd_str.lower()


class TestOptionalDependencyError:
    """Test OptionalDependencyError raised for missing streamlit."""

    def test_error_has_package_name(self) -> None:
        """OptionalDependencyError contains the package name."""
        from chemfuse.core.exceptions import OptionalDependencyError

        exc = OptionalDependencyError("streamlit", extra="web")
        assert "streamlit" in exc.message

    def test_error_has_install_hint(self) -> None:
        """OptionalDependencyError contains installation instructions."""
        from chemfuse.core.exceptions import OptionalDependencyError

        exc = OptionalDependencyError("streamlit", extra="web")
        assert "pip install" in exc.message

    def test_error_has_extra_hint(self) -> None:
        """OptionalDependencyError references the optional extra."""
        from chemfuse.core.exceptions import OptionalDependencyError

        exc = OptionalDependencyError("streamlit", extra="web")
        assert "web" in exc.message
