"""ChemFuse CLI main entry point."""

from __future__ import annotations

import click

from chemfuse._version import __version__
from chemfuse.cli.commands.profile import profile_cmd
from chemfuse.cli.commands.screen import screen_cmd
from chemfuse.cli.commands.search import search_cmd
from chemfuse.cli.commands.serve import serve
from chemfuse.cli.commands.web import web_cmd
from chemfuse.cli.commands.xref import xref_cmd


@click.group()
@click.version_option(version=__version__, prog_name="chemfuse")
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help="Enable verbose output.",
)
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """ChemFuse: Multi-database cheminformatics suite.

    Search and analyze compounds across PubChem and other chemical databases.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


cli.add_command(search_cmd)
cli.add_command(screen_cmd)
cli.add_command(profile_cmd)
cli.add_command(xref_cmd)
cli.add_command(web_cmd)
cli.add_command(serve)


if __name__ == "__main__":
    cli()
