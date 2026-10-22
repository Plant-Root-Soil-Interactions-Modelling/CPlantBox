"""Erzeugt aus der Template-XML ein Parameterfile für einen Designpunkt."""
import xml.etree.ElementTree as ET
from pathlib import Path

import config


def _parse(path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))  # Kommentare behalten
    return ET.parse(path, parser)


def _hyphae(root, subtype):
    for el in root.iter("hyphae"):
        if el.get("subType", "").strip() == str(subtype):
            return el
    raise KeyError(f"Kein <hyphae subType='{subtype}'> in der Template-XML")


def _param(hy, name):
    for p in hy.iter("parameter"):
        if p.get("name") == name:
            return p
    raise KeyError(f"Parameter '{name}' fehlt in subType {hy.get('subType')}")


def _set_value(p, new_value):
    old = float(p.get("value"))
    p.set("value", repr(float(new_value)))
    if p.get("dev") is not None and config.DEV_MODE == "relative" and old != 0:
        p.set("dev", repr(float(p.get("dev")) * new_value / old))


def _set_successor_p(hy, p_first):
    succ = next(q for q in hy.iter("parameter") if q.get("name") == "successor")
    succ.set("percentage", f"{p_first:.6g}, {1 - p_first:.6g}")


def write_xml(values: dict, out_path: Path, dx: float | None = None, template=None):
    """values: {factor_name: physical value}. Nicht genannte Faktoren bleiben nominal."""
    tree = _parse(template or config.TEMPLATE_XML)
    root = tree.getroot()
    fmap = {f.name: f for f in config.FACTORS}
    for name, val in values.items():
        f = fmap[name]
        hy = _hyphae(root, f.subtype)
        if f.param == "percentage":
            _set_successor_p(hy, val)
        else:
            _set_value(_param(hy, f.param), val)
    if dx is not None:                       # nur Phase 0 (Konvergenztest)
        for hy in root.iter("hyphae"):
            _param(hy, "dx").set("value", repr(float(dx)))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    return out_path


def nominal_values(template=None):
    """Liest die Nominalwerte aller Faktoren aus der Template-XML."""
    root = _parse(template or config.TEMPLATE_XML).getroot()
    out = {}
    for f in config.FACTORS:
        hy = _hyphae(root, f.subtype)
        if f.param == "percentage":
            succ = next(q for q in hy.iter("parameter") if q.get("name") == "successor")
            out[f.name] = float(succ.get("percentage").split(",")[0])
        else:
            out[f.name] = float(_param(hy, f.param).get("value"))
    return out
