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


def _set_cv(p, cv):
    p.set("dev", repr(float(cv) * float(p.get("value"))))


def write_xml(values: dict, out_path: Path, dx: float | None = None, template=None):
    """values: {factor_name: physical value}. Nicht genannte Faktoren bleiben nominal."""
    tree = _parse(template or config.TEMPLATE_XML)
    root = tree.getroot()
    fmap = {f.name: f for f in config.FACTORS}
    # erst alle Werte, dann v_cv: dev von v bezieht sich auf das bereits gesetzte v
    for name in sorted(values, key=lambda n: fmap[n].param == "v_cv"):
        f, val = fmap[name], values[name]
        for st in f.subtypes:
            hy = _hyphae(root, st)
            if f.param == "v_cv":
                _set_cv(_param(hy, "v"), val)
            else:
                _set_value(_param(hy, f.param), val)
    if dx is not None:                       # nur Phase 0 (Konvergenztest)
        for hy in root.iter("hyphae"):
            _param(hy, "dx").set("value", repr(float(dx)))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    return out_path


def nominal_values(template=None):
    """Liest die Nominalwerte aller Faktoren aus der Template-XML (gemeinsame Faktoren: erster subType)."""
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
