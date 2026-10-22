"""
Prüft die Konfiguration, bevor irgendetwas auf den Cluster geht.

    python check_config.py
"""
import xml.etree.ElementTree as ET

import config
from xml_tools import nominal_values

MIN_PER_RUN = 20  # Simulation + Auswertung [min], für die Abschätzung


def main():
    nom = nominal_values()
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
        if f.param in ("ana", "percentage") and not (0 <= f.low and f.high <= 1):
            problems.append(f"{f.name}: Wahrscheinlichkeit muss in [0, 1] liegen")

    # b_prob muss 0/1 sein, damit die Faktorauswahl stimmt
    root = ET.parse(config.TEMPLATE_XML).getroot()
    for hy in root.iter("hyphae"):
        bp = next((p.get("value") for p in hy.iter("parameter") if p.get("name") == "b_prob"), None)
        print(f"\nsubType {hy.get('subType')} ({hy.get('name')}): b_prob = {bp}", end="")
    print()

    k = len(config.FACTORS)
    n0 = config.PHASE0_NOISE_SEEDS + len(config.PHASE0_DX_VALUES) * config.PHASE0_DX_SEEDS
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
