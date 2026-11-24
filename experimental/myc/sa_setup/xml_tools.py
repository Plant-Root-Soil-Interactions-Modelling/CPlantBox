"""Creates a parameter file for a design point from the template XML."""
import xml.etree.ElementTree as ET
from pathlib import Path

import config


def _parse(path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))  # keep comments
    return ET.parse(path, parser)


def _hyphae(root, subtype):
    for el in root.iter("hyphae"):
        if el.get("subType", "").strip() == str(subtype):
            return el
    raise KeyError(f"No <hyphae subType='{subtype}'> in the template XML")


def _param(hy, name):
    for p in hy.iter("parameter"):
        if p.get("name") == name:
            return p
    raise KeyError(f"Parameter '{name}' missing in subType {hy.get('subType')}")


def _set_value(p, new_value):
    old = float(p.get("value"))
    p.set("value", repr(float(new_value)))
    if p.get("dev") is not None and config.DEV_MODE == "relative" and old != 0:
        p.set("dev", repr(float(p.get("dev")) * new_value / old))


def _set_cv(p, cv):
    p.set("dev", repr(float(cv) * float(p.get("value"))))


def write_xml(values: dict, out_path: Path, dx: float | None = None, template=None):
    """values: {factor_name: physical value}. Factors not given stay nominal."""
    tree = _parse(template or config.TEMPLATE_XML)
    root = tree.getroot()
    fmap = {f.name: f for f in config.FACTORS}
    # all values first, then v_cv: dev of v refers to the v already set
    for name in sorted(values, key=lambda n: fmap[n].param == "v_cv"):
        f, val = fmap[name], values[name]
        for st in f.subtypes:
            hy = _hyphae(root, st)
            if f.param == "v_cv":
                _set_cv(_param(hy, "v"), val)
            else:
                _set_value(_param(hy, f.param), val)
    if dx is not None:                       # Phase 0 only (convergence test)
        for hy in root.iter("hyphae"):
            _param(hy, "dx").set("value", repr(float(dx)))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    return out_path


def write_hyphae_xml(updates: dict, out_path: Path, template=None):
    """updates: {(subType, parameter name): value}, independent of config.FACTORS (calibration).
    "successor" sets the share of the first successor: percentage = "p, 1 - p"."""
    tree = _parse(template or config.TEMPLATE_XML)
    root = tree.getroot()
    for (st, name), val in updates.items():
        p = _param(_hyphae(root, st), name)
        if name == "successor":
            if len(p.get("percentage").split(",")) != 2:
                raise ValueError(f"successor of subType {st}: exactly two successors expected")
            p.set("percentage", f"{float(val)!r}, {1 - float(val)!r}")
        else:
            _set_value(p, val)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    return out_path


def hyphae_nominal(subtype, name, template=None):
    """Nominal value of a hyphal parameter (successor: share of the first successor)."""
    p = _param(_hyphae(_parse(template or config.TEMPLATE_XML).getroot(), subtype), name)
    return float(p.get("percentage").split(",")[0]) if name == "successor" else float(p.get("value"))


def nominal_values(template=None):
    """Reads the nominal values of all factors from the template XML (shared factors: first subType)."""
    root = _parse(template or config.TEMPLATE_XML).getroot()
    out = {}
    for f in config.FACTORS:
        hy = _hyphae(root, f.subtypes[0])
        if f.param == "v_cv":
            p = _param(hy, "v")
            out[f.name] = float(p.get("dev") or 0) / float(p.get("value"))
        else:
            out[f.name] = float(_param(hy, f.param).get("value"))
    return out
