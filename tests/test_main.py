"""Tests for the __main__ module entry point."""


def test_main_module_importable() -> None:
    """Verify that the __main__ entry point can be imported without error."""
    import chemfuse.__main__  # noqa: F401
