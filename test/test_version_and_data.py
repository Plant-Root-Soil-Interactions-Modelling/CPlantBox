import os
import re


def test_version_pep440():
    import plantbox as pb

    v = getattr(pb, "__version__", "0+unknown")
    # Allow unknown in dev trees; otherwise basic PEP 440-compatible shape
    if v != "0+unknown":
        assert re.match(r"^\d+\.\d+(\.\d+)?([ab]|rc)?\d*(\.post\d+)?(\.dev\d+)?(\+.*)?$", v)


def test_data_path_exists():
    import plantbox as pb

    dp = pb.data_path()
    assert os.path.isdir(dp)
    # A structural subdir should exist
    assert os.path.isdir(os.path.join(dp, "structural"))
