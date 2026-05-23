#!/usr/bin/env python3
"""Smoke-test an installed CPlantBox wheel from an end-user perspective."""

from __future__ import annotations

from pathlib import Path
import importlib.util
import site
import sys
import sysconfig


REPO_ROOT = Path(__file__).resolve().parents[2]


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _site_package_paths() -> list[Path]:
    paths: list[Path] = []
    for key in ("platlib", "purelib"):
        value = sysconfig.get_paths().get(key)
        if value:
            paths.append(Path(value).resolve())
    try:
        paths.extend(Path(value).resolve() for value in site.getsitepackages())
    except AttributeError:
        pass
    return list(dict.fromkeys(paths))


def _prefer_installed_package() -> list[Path]:
    """Put site-packages first and remove CPlantBox source-tree paths."""

    site_paths = _site_package_paths()
    sanitized: list[str] = []
    cwd = Path.cwd().resolve()

    for entry in sys.path:
        resolved = cwd if entry == "" else Path(entry).resolve()
        if resolved == REPO_ROOT or _is_relative_to(resolved, REPO_ROOT):
            continue
        sanitized.append(entry)

    site_strings = [str(path) for path in site_paths]
    sys.path[:] = site_strings + [entry for entry in sanitized if entry not in site_strings]
    return site_paths


def _assert_installed_path(path: Path, site_paths: list[Path]) -> None:
    resolved = path.resolve()
    assert not _is_relative_to(resolved, REPO_ROOT), resolved
    if site_paths:
        assert any(_is_relative_to(resolved, site_path) for site_path in site_paths), (resolved, site_paths)


def _assert_module_packaged(module_name: str, site_paths: list[Path]) -> None:
    spec = importlib.util.find_spec(module_name)
    assert spec is not None, module_name
    assert spec.origin is not None, module_name
    _assert_installed_path(Path(spec.origin), site_paths)


def _check_helper_modules(site_paths: list[Path]) -> None:
    """Check representative installed helper modules and subpackages."""

    # Import helpers covered by the base runtime deps plus cibuildwheel's `vis` test extra.
    import plantbox.functional.Perirhizal  # noqa: F401
    import plantbox.functional.PlantHydraulicModel  # noqa: F401
    import plantbox.functional.van_genuchten  # noqa: F401
    import plantbox.rsml.rsml_reader  # noqa: F401
    import plantbox.rsml.rsml_writer  # noqa: F401
    import plantbox.structural.MappedOrganism  # noqa: F401
    import plantbox.visualisation.vtk_plot  # noqa: F401
    import plantbox.visualisation.vtk_tools  # noqa: F401

    # FiPy helper modules import mpi4py, which requires an external MPI runtime
    # library. Verify they are packaged without executing those optional imports.
    _assert_module_packaged("plantbox.functional.CellVariablemod", site_paths)
    _assert_module_packaged("plantbox.functional.Mesh1Dmod", site_paths)


def main() -> None:
    site_paths = _prefer_installed_package()

    import plantbox as pb
    import plantbox.functional.PerirhizalHeterogeneous as perirhizal_heterogeneous
    import plantbox.plantbox as extension

    assert hasattr(perirhizal_heterogeneous, "PerirhizalHetereogeneous")

    module_path = Path(pb.__file__).resolve()
    extension_path = Path(extension.__file__).resolve()
    _assert_installed_path(module_path, site_paths)
    _assert_installed_path(extension_path, site_paths)

    _check_helper_modules(site_paths)

    data_root = Path(pb.data_path()).resolve()
    _assert_installed_path(data_root, site_paths)
    assert data_root.is_dir(), data_root

    parameter_file = data_root / "structural" / "plant" / "fspm2023.xml"
    assert parameter_file.is_file(), parameter_file

    plant = pb.Plant(1)
    plant.readParameters(str(parameter_file))
    plant.initialize(False)
    plant.simulate(1, False)

    assert plant.getNumberOfOrgans() > 0
    assert plant.getNumberOfNodes() > 0
    assert plant.getNumberOfSegments() > 0

    print(f"plantbox package: {module_path}")
    print(f"plantbox extension: {extension_path}")
    print(f"model parameters: {data_root}")
    print(
        "smoke result: "
        f"{plant.getNumberOfOrgans()} organs, "
        f"{plant.getNumberOfNodes()} nodes, "
        f"{plant.getNumberOfSegments()} segments"
    )


if __name__ == "__main__":
    main()
