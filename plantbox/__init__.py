from importlib import import_module as _import_module

# Try to import the binary submodule shipped inside the package dir first
try:
    _ext = _import_module("._plantbox", package=__name__)
except Exception as _e:
    try:
        _ext = _import_module("_plantbox")
    except Exception:  # pragma: no cover
        raise

for _name in dir(_ext):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_ext, _name)

__all__ = [n for n in dir(_ext) if not n.startswith("_")]

# Expose package version from installed metadata (falls back to bundled _version.py)
try:  # pragma: no cover
    # Prefer the wheel metadata version to avoid SCM/runtime dependencies
    from importlib.metadata import version as _pkg_version  # type: ignore

    __version__ = _pkg_version("cplantbox")  # project name
except Exception:  # pragma: no cover
    try:
        from ._version import version as __version__  # type: ignore
    except Exception:
        __version__ = "0+unknown"


def data_path() -> str:
    """Return the root directory for packaged model parameter data.

    Works from both an installed wheel (where data is under the package dir)
    and a source checkout (falls back to repository tree).
    """
    import os

    pkg_dir = os.path.dirname(__file__)
    # First preference: packaged under the module directory (wheel install)
    candidate = os.path.join(pkg_dir, "modelparameter")
    if os.path.isdir(candidate):
        return candidate
    # Fallback: source tree (repo layout)
    # Two levels up from package dir should be repo root if run in-place
    repo_root = os.path.abspath(os.path.join(pkg_dir, os.pardir))
    fallback = os.path.join(repo_root, "modelparameter")
    return fallback
