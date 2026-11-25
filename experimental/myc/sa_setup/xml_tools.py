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


def _scale_roots(root, name, factor):
    """Multiplies the <parameter> `name` (value and dev) of all <root> elements by factor."""
    n = 0
    for el in root.iter("root"):
        for p in el.iter("parameter"):
            if p.get("name") == name:
                p.set("value", repr(float(p.get("value")) * factor))
                if p.get("dev") is not None:
                    p.set("dev", repr(float(p.get("dev")) * factor))
                n += 1
    if n == 0:
        raise KeyError(f"No <root> parameter '{name}' in the template XML")


def derived_ana_min_dist(hy):
    """Minimal path length for anastomosis as the model derives it for anaMinDist < 0: distTH / sin(theta) + 2 dx [cm]."""
    import math
    v = lambda n: float(_param(hy, n).get("value"))
    return v("distTH") / max(math.sin(v("theta")), 1e-6) + 2 * v("dx")


def _set_ana_L(hy, scale):
    """anaMinDist = scale × derived value (from distTH, theta, dx of this <hyphae> as already written)."""
    try:
        p = _param(hy, "anaMinDist")
    except KeyError:
        p = ET.SubElement(hy, "parameter", name="anaMinDist")
    p.set("value", repr(float(scale) * derived_ana_min_dist(hy)))


def _set_successor(hy, share):
    """Share of the first successor: percentage = "p, 1 - p"."""
    p = _param(hy, "successor")
    if len(p.get("percentage").split(",")) != 2:
        raise ValueError(f"successor of subType {hy.get('subType')}: exactly two successors expected")
    p.set("percentage", f"{float(share)!r}, {1 - float(share)!r}")


def _scale_devs(hy, scale):
    """All dev of this <hyphae> except the one of v (v_cv) × scale (cv_scale, A9)."""
    for p in hy.iter("parameter"):
        if p.get("name") != "v" and p.get("dev") is not None:
            p.set("dev", repr(float(p.get("dev")) * float(scale)))


def write_xml(values: dict, out_path: Path, dx: float | None = None, template=None, ana_L_scale: float | None = None,
              factors=None):
    """values: {factor_name: physical value}. factors: their definitions (default config.FACTORS + INOCULUM_FACTORS; the
    first Morris round: config.FACTORS_V1). Factors not given stay nominal. Factors with target "sim" are not written
    (they go to model_hook.simulate). Special params: "successor" (share of the first successor), "v_cv" (dev of v),
    "ana_L_scale" (anaMinDist = scale × (distTH / sin(theta) + 2 dx), after distTH, theta; without it the template's
    anaMinDist = -1 lets the model derive it), "cv_scale" (all dev but v, last). `also`: further parameters with the same value.
    ana_L_scale argument: Phase 0 test of the minimal path length."""
    tree = _parse(template or config.TEMPLATE_XML)
    root = tree.getroot()
    fmap = {f.name: f for f in (factors if factors is not None else config.FACTORS) + config.INOCULUM_FACTORS}
    order = {"v_cv": 1, "ana_L_scale": 2, "cv_scale": 3}   # dev of v refers to v; anaMinDist to distTH, theta
    for name in sorted(values, key=lambda n: order.get(fmap[n].param, 0)):
        f, val = fmap[name], values[name]
        if f.target == "root_scale":
            _scale_roots(root, f.param, val)
            continue
        if f.target == "sim":
            continue
        for st in f.subtypes:
            hy = _hyphae(root, st)
            if f.param == "v_cv":
                _set_cv(_param(hy, "v"), val)
            elif f.param == "ana_L_scale":
                _set_ana_L(hy, val)
            elif f.param == "cv_scale":
                _scale_devs(hy, val)
            elif f.param == "successor":
                _set_successor(hy, val)
            else:
                for pname in (f.param,) + tuple(f.also):
                    _set_value(_param(hy, pname), val)
    if dx is not None:                       # Phase 0 only (convergence test)
        for hy in root.iter("hyphae"):
            _param(hy, "dx").set("value", repr(float(dx)))
    if ana_L_scale is not None:              # Phase 0 only (test of the minimal path length)
        for hy in root.iter("hyphae"):
            _set_ana_L(hy, ana_L_scale)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    return out_path


def hyphae_nominal(subtype, name, template=None):
    """Nominal value of a hyphal parameter (successor: share of the first successor)."""
    p = _param(_hyphae(_parse(template or config.TEMPLATE_XML).getroot(), subtype), name)
    return float(p.get("percentage").split(",")[0]) if name == "successor" else float(p.get("value"))


def factor_nominal(f, template=None):
    """Nominal value of a factor: Factor.nominal if given, v_cv = dev/value of v, else from the template XML (first subType)."""
    if f.nominal is not None:
        return float(f.nominal)
    if f.param == "v_cv":
        p = _param(_hyphae(_parse(template or config.TEMPLATE_XML).getroot(), f.subtypes[0]), "v")
        return float(p.get("dev") or 0) / float(p.get("value"))
    return hyphae_nominal(f.subtypes[0], f.param, template)


def nominal_values(template=None, factors=None):
    """Nominal values of all factors (default config.FACTORS; shared factors: first subType)."""
    return {f.name: factor_nominal(f, template) for f in (factors if factors is not None else config.FACTORS)}
