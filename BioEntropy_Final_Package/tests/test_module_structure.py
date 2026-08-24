"""Architecture guards for the module split.

These pin the dependency layering so the extracted modules don't silently grow
a back-dependency on the heavy core (which would reintroduce import cycles).
"""
import ast
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent


def _imported_modules(filename: str) -> set[str]:
    tree = ast.parse((PKG / filename).read_text(encoding="utf-8"))
    mods: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module)
    return mods


def test_constants_is_a_leaf():
    # The constants module must not import any project module.
    project = {m for m in _imported_modules("bioentropy_constants.py")
               if m.startswith("bioentropy")}
    assert project == set(), f"constants should be dependency-free, got {project}"


def test_entropy_only_depends_on_constants():
    # The numerical primitives must not import the heavy core (would be a cycle).
    project = {m for m in _imported_modules("bioentropy_entropy.py")
               if m.startswith("bioentropy")}
    assert project == {"bioentropy_constants"}, (
        f"entropy must depend only on bioentropy_constants, got {project}"
    )


def test_display_state_is_a_leaf():
    project = {m for m in _imported_modules("bioentropy_display_state.py")
               if m.startswith("bioentropy")}
    assert project == set(), f"display_state should be dependency-free, got {project}"


def test_plots_does_not_import_the_orchestrator():
    # bioentropy_core_generalized imports bioentropy_plots; the reverse would be
    # a cycle. Plots may depend on the legacy core, constants and display state.
    project = {m for m in _imported_modules("bioentropy_plots.py")
               if m.startswith("bioentropy")}
    assert "bioentropy_core_generalized" not in project, (
        f"plots must not import the orchestrator, got {project}"
    )
    assert project <= {"bioentropy_core_legacy", "bioentropy_constants",
                       "bioentropy_display_state"}, f"unexpected plots deps: {project}"


def test_modules_import_standalone():
    import bioentropy_constants  # noqa: F401
    import bioentropy_entropy  # noqa: F401
    import bioentropy_display_state  # noqa: F401
    import bioentropy_plots  # noqa: F401

    assert bioentropy_constants.SOFTWARE_VERSION
    # numerical primitives are importable without the heavy core
    assert hasattr(bioentropy_entropy, "compute_group_entropy")
    assert hasattr(bioentropy_entropy, "infer_direction")
    # plotting entry points
    assert hasattr(bioentropy_plots, "plot_entropy_trajectory")
