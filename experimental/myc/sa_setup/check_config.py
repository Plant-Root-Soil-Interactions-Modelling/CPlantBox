"""
Prüft die Konfiguration, bevor irgendetwas auf den Cluster geht.

    python check_config.py
"""
import xml.etree.ElementTree as ET

import config
from xml_tools import nominal_values


def nominal_values_of(param, subtype):
    root = ET.parse(config.TEMPLATE_XML).getroot()
    hy = next(h for h in root.iter("hyphae") if h.get("subType", "").strip() == str(subtype))
    return float(next(p for p in hy.iter("parameter") if p.get("name") == param).get("value"))

MIN_PER_RUN = 120  # Simulation + Auswertung [min], für die Abschätzung (Seeds brauchen bis zu 3 h, s. petridishtiming.py)


def main():
    nom = nominal_values()
    print(f"Nominalmodell: {config.TEMPLATE_XML}")
    print(f"{'Faktor':10s} {'nominal':>10s} {'low':>10s} {'high':>10s}  Quelle")
    open_src, problems = 0, []
    for f in config.FACTORS:
        n = nom[f.name]
        src = f.source or "OFFEN (Platzhalter)"
        open_src += f.source is None
        print(f"{f.name:10s} {n:10.4g} {f.low:10.4g} {f.high:10.4g}  {src}"
              + ("   [Kontrolle]" if f.role == "control" else ""))
        if not (f.low <= n <= f.high):
            problems.append(f"{f.name}: Nominalwert {n} liegt außerhalb [{f.low}, {f.high}]")
        if f.low >= f.high:
            problems.append(f"{f.name}: low >= high")
        if f.log and f.low <= 0:
            problems.append(f"{f.name}: log-Bereich braucht low > 0")
        if f.param == "ana" and not (0 <= f.low and f.high <= 1):
            problems.append(f"{f.name}: Wahrscheinlichkeit muss in [0, 1] liegen")
    # gemeinsame Faktoren: alle subTypes sollten denselben Nominalwert haben
    for f in config.FACTORS:
        if len(f.subtypes) > 1 and f.param != "v_cv":
            vals = {st: nominal_values_of(f.param, st) for st in f.subtypes}
            if len(set(vals.values())) > 1:
                problems.append(f"{f.name}: gemeinsamer Faktor, aber verschiedene Nominalwerte {vals}")

    # b_prob muss 0/1 sein, damit die Faktorauswahl stimmt
    root = ET.parse(config.TEMPLATE_XML).getroot()
    for hy in root.iter("hyphae"):
        bp = next((p.get("value") for p in hy.iter("parameter") if p.get("name") == "b_prob"), None)
        print(f"\nsubType {hy.get('subType')} ({hy.get('name')}): b_prob = {bp}", end="")
    print()

    k = len(config.FACTORS)
    n0 = config.PHASE0_NOISE_SEEDS + len(config.PHASE0_DX_VALUES) * config.PHASE0_DX_SEEDS + config.PHASE0_RANGE_POINTS
    n1 = config.MORRIS_TRAJECTORIES * (k + 1) * config.MORRIS_REPLICATES
    print(f"\nPhase 0: {n0} Läufe  ≈ {n0 * MIN_PER_RUN / 60:.0f} CPU-h")
    print(f"Morris:  {config.MORRIS_TRAJECTORIES} Trajektorien × {k + 1} Punkte × "
          f"{config.MORRIS_REPLICATES} Replikate = {n1} Läufe  ≈ {n1 * MIN_PER_RUN / 60:.0f} CPU-h")
    print(f"dev-Modus: {config.DEV_MODE}")
    if open_src:
        print(f"\n{open_src} von {k} Faktorbereichen sind noch Platzhalter (±25 %).")
    for p in problems:
        print("FEHLER:", p)


if __name__ == "__main__":
    main()
