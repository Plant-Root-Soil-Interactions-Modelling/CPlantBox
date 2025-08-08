import os
import shutil

import pytest


def _has_display() -> bool:
    return bool(
        os.environ.get("DISPLAY")
        or os.environ.get("WAYLAND_DISPLAY")
        or os.environ.get("XDG_SESSION_TYPE") == "wayland"
    )


def _has_vtk() -> bool:
    try:
        import vtk  # noqa: F401

        return True
    except Exception:
        return False


def pytest_collection_modifyitems(config, items):
    vtk_available = _has_vtk()
    display_available = _has_display()
    skip_vtk = pytest.mark.skip(reason="requires VTK which is not available in this environment")
    skip_gui = pytest.mark.skip(
        reason="requires GUI/display which is not available in this environment"
    )
    for item in items:
        if "vtk" in item.keywords and not vtk_available:
            item.add_marker(skip_vtk)
        if "gui" in item.keywords and not display_available:
            item.add_marker(skip_gui)


@pytest.fixture(autouse=True)
def _chdir_tmp_and_cleanup(tmp_path, request):
    """Run each test in its own temp directory and clean legacy artifacts.

    Many legacy tests write fixed filenames in CWD. This fixture changes CWD
    to a unique per-test temp dir to avoid polluting the repo. As a fallback
    (for tests that still insist on writing into repo root) we remove known
    legacy artifacts post-test if they were created there.
    """
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(old_cwd)

    # Fallback cleanup in repo root for known legacy outputs
    legacy = [
        "human.xml",
        "leaf.xml",
        "organ.xml",
        "organism.rsml",
        "root.xml",
        "rs_parameters.xml",
        "seed.xml",
        "stem.xml",
        "test_rootsystem.rsml",
        "test_rotsystem.vtp",
        "test_rootsystem.vtp",
    ]
    for name in legacy:
        path = os.path.join(repo_root, name)
        try:
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)
        except Exception:
            # Best-effort cleanup; do not fail tests on cleanup errors
            pass
