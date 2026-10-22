"""
Zentrale Konfiguration der Sensitivitätsanalyse.

Hier wird festgelegt
  - welche Parameter (Faktoren) variiert werden und in welchem Bereich,
  - wie mit den `dev`-Werten umgegangen wird,
  - Seeds, Replikate, Zeitpunkte, Zielbereich und Teilbereiche.

ALLE BEREICHE SIND PLATZHALTER (±25 % um den Nominalwert), solange keine
Literatur-/Datenbasis eingetragen ist. Faktoren mit `source=None` werden in
`check_config.py` als offen gemeldet.
"""
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Pfade
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
TEMPLATE_XML = ROOT / "template" / "hyphae_template.xml"   # dein Nominal-Parameterfile
RUNS_DIR = ROOT / "runs"          # pro Lauf ein Ordner mit XML, Rohausgabe, meta.json
RESULTS_DIR = ROOT / "results"    # Designs, Output-Tabellen, Analysen, Plots

# --------------------------------------------------------------------------
# Umgang mit `dev`
# --------------------------------------------------------------------------
# "relative": dev wird mit dem Mittelwert skaliert (CV = dev/value bleibt konstant).
#             Verhindert z. B., dass lmax=3.5 mit dev=3 massenhaft negative Ziehungen erzeugt.
# "absolute": dev bleibt auf dem Nominalwert aus der Template-XML.
DEV_MODE = "relative"

# --------------------------------------------------------------------------
# Faktoren
# --------------------------------------------------------------------------
@dataclass
class Factor:
    name: str            # eindeutiger Name in Design/Analyse, z. B. "P_v"
    subtype: int         # subType-Attribut des <hyphae>-Elements
    param: str           # name-Attribut des <parameter>-Elements (oder "percentage")
    low: float
    high: float
    log: bool = False    # True -> log-uniform zwischen low und high
    source: str | None = None   # Begründung für den Bereich (Literatur, Daten, Annahme)
    role: str = "factor"        # "factor" oder "control" (sollte ~0 Einfluss haben)


def _pm(rel, nominal):
    """Platzhalterbereich nominal*(1±rel)."""
    return nominal * (1 - rel), nominal * (1 + rel)


PH = 0.25  # Platzhalter-Spannweite ±25 %

PULLER, BAS = 1, 2

FACTORS: list[Factor] = [
    # ---- PullerHyphae (subType 1, laterale Verzweigung, b_prob = 0) ----
    Factor("P_a",      PULLER, "a",      *_pm(PH, 0.005)),
    Factor("P_ana",    PULLER, "ana",    0.5, 1.0),           # Wahrscheinlichkeit, nominal 1 = Rand
    Factor("P_distTH", PULLER, "distTH", *_pm(PH, 0.02)),
    Factor("P_hlt",    PULLER, "hlt",    *_pm(PH, 10)),
    Factor("P_la",     PULLER, "la",     *_pm(PH, 0.01)),
    Factor("P_lb",     PULLER, "lb",     *_pm(PH, 0.4)),
    Factor("P_lmax",   PULLER, "lmax",   *_pm(PH, 7)),
    Factor("P_ln",     PULLER, "ln",     *_pm(PH, 0.1)),
    Factor("P_theta",  PULLER, "theta",  *_pm(PH, 1.0471976)),
    Factor("P_v",      PULLER, "v",      *_pm(PH, 0.68)),
    # Kontrollfaktor: b ist bei b_prob = 0 wirkungslos -> erwartet: mu* ~ Rauschniveau
    Factor("P_b_ctrl", PULLER, "b",      *_pm(PH, 0.95), role="control"),

    # ---- BAS (subType 2, Tip-Splitting, b_prob = 1) ----
    Factor("BAS_a",      BAS, "a",      *_pm(PH, 0.005)),
    Factor("BAS_ana",    BAS, "ana",    0.5, 1.0),
    Factor("BAS_b",      BAS, "b",      *_pm(PH, 24)),
    Factor("BAS_distTH", BAS, "distTH", *_pm(PH, 0.02)),
    Factor("BAS_hlt",    BAS, "hlt",    *_pm(PH, 6)),
    Factor("BAS_lmax",   BAS, "lmax",   *_pm(PH, 3)),
    Factor("BAS_theta",  BAS, "theta",  *_pm(PH, 1.0471976)),
    Factor("BAS_v",      BAS, "v",      *_pm(PH, 0.60)),

    # ---- Successor-Aufteilung der Puller-Laterals: "p, 1-p" für subType "2, 1" ----
    Factor("P_succ_p", PULLER, "percentage", 0.5, 1.0),       # nominal 0.75
]

# Parameter, die NIE variiert werden (nur zur Dokumentation / Prüfung)
FIXED = ["tropismT", "tropismN", "tropismS", "dxMin", "ldelay", "b_prob", "dx"]

# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------
SIM_TIME = 30                     # Simulationsdauer [Tage]  TODO: anpassen
SNAPSHOT_TIMES = [10, 20, 30]     # Zeitpunkte, zu denen Rohausgabe gespeichert wird

# Zielbereich (Box) und Unterteilung. Achsen: x, y, z [cm]   TODO: anpassen
DOMAIN = {"min": (-5.0, -5.0, -10.0), "max": (5.0, 5.0, 0.0)}
REGIONS = {
    # Teilbereiche als Boxen. Default: 5 Tiefenschichten à 2 cm.
    f"z{i}": {"min": (-5.0, -5.0, -2.0 * (i + 1)), "max": (5.0, 5.0, -2.0 * i)}
    for i in range(5)
}

# --------------------------------------------------------------------------
# Phase 0
# --------------------------------------------------------------------------
PHASE0_NOISE_SEEDS = 50
PHASE0_DX_VALUES = [0.01, 0.005, 0.0025]   # nominal, halb, viertel
PHASE0_DX_SEEDS = 10

# --------------------------------------------------------------------------
# Phase 1 (Morris)
# --------------------------------------------------------------------------
MORRIS_LEVELS = 4            # p
MORRIS_TRAJECTORIES = 12     # r (aus MORRIS_CANDIDATES nach Spreizung ausgewählt)
MORRIS_CANDIDATES = 300
MORRIS_REPLICATES = 3        # nach Phase 0 anpassen
BASE_SEED = 20261021
