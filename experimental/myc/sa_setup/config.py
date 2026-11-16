"""
Zentrale Konfiguration der Sensitivitätsanalyse (Petrischalen-Setup, TwoHyphaePlusBAS.xml).

Hier wird festgelegt
  - welche Parameter (Faktoren) variiert werden und in welchem Bereich,
  - wie mit den `dev`-Werten umgegangen wird,
  - Seeds, Replikate, Simulationsdauer nach t*, Auswertezeitpunkte und Outputs.

Bereiche mit `source=None` sind Platzhalter (±25 % um den Nominalwert) und werden in
`check_config.py` als offen gemeldet.
"""
import math
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Pfade
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
MYC_DIR = ROOT.parent                              # experimental/myc (AMF*.py)
CPLANTBOX_DIR = ROOT.parents[2]                    # für `import plantbox`
TEMPLATE_XML = MYC_DIR / "tomatoparameters" / "TwoHyphaePlusBAS.xml"   # Nominalmodell
RUNS_DIR = ROOT / "runs"          # pro Lauf ein Ordner mit XML, Rohausgabe, meta.json
RESULTS_DIR = ROOT / "results"    # Designs, Output-Tabellen, Analysen, Plots

# --------------------------------------------------------------------------
# Umgang mit `dev`
# --------------------------------------------------------------------------
# "relative": dev wird mit dem Mittelwert skaliert (CV = dev/value bleibt konstant).
# "absolute": dev bleibt auf dem Nominalwert aus der Template-XML.
# Der Faktor `v_cv` setzt dev von v unabhängig davon auf CV × v.
DEV_MODE = "relative"

# --------------------------------------------------------------------------
# Faktoren
# --------------------------------------------------------------------------
RUNNER, BAS = 1, 2   # subTypes in TwoHyphaePlusBAS.xml
SUBTYPE_NAMES = {RUNNER: "Runner", BAS: "BAS"}


@dataclass
class Factor:
    name: str               # eindeutiger Name in Design/Analyse, z. B. "R_v"
    subtypes: tuple         # subTypes der <hyphae>, auf die der Wert geschrieben wird (mehrere = gemeinsamer Faktor)
    param: str              # name-Attribut des <parameter>-Elements, oder "v_cv" (dev von v = CV × v)
    low: float
    high: float
    log: bool = False       # True -> log-uniform zwischen low und high
    source: str | None = None   # Begründung für den Bereich (Literatur, Daten, Annahme)
    role: str = "factor"        # "factor" oder "control" (sollte keinen Einfluss haben)


def _pm(rel, nominal):
    """Platzhalterbereich nominal*(1±rel)."""
    return nominal * (1 - rel), nominal * (1 + rel)


PH = 0.25  # Platzhalter-Spannweite ±25 %
V_SOURCE = "240 µm/h (Oyarte Galvez et al. 2025), Bereich 180–300 µm/h um den Schwerpunkt der gemessenen Spitzengeschwindigkeiten"

FACTORS: list[Factor] = [
    # ---- RunnerHyphae (subType 1, laterale Verzweigung, b_prob = 0) ----
    Factor("R_v",      (RUNNER,), "v",     0.432, 0.72, source=V_SOURCE),
    Factor("R_hlt",    (RUNNER,), "hlt",   4.0, 10.0,
           source="Sensitivitätsbereich Schnepf 2016; Sterberaten Schnepf 2008 ergeben 4–12,5 d"),
    Factor("R_lmax",   (RUNNER,), "lmax",  1.5, 5.0,
           source="nicht in der Literatur; Abschätzung v × Lebensdauer (≈ 1,3–5,4 cm)"),
    Factor("R_ln",     (RUNNER,), "ln",    0.1, 0.6, log=True,
           source="~1 mm alle Äste (Giovannetti, in Schnepf 2008) bis ~6 mm neue wachsende Spitzen "
                  "(Oyarte Galvez 2025); log-uniform, da Faktor 6 und Lateralrate ∝ 1/ln"),
    Factor("R_theta",  (RUNNER,), "theta", math.radians(50), math.radians(80),
           source="60° (Schnepf 2016) bis ≈ 70° (Oyarte Galvez 2025), Winkeldefinitionen unterscheiden sich"),
    # Kontrollfaktor: b wird nur beim Tip-Splitting gelesen (Hyphae.cpp), bei b_prob = 0 wirkungslos -> mu* = 0
    Factor("R_b_ctrl", (RUNNER,), "b",     *_pm(PH, 0.95), role="control"),

    # ---- BAS (subType 2, Tip-Splitting, b_prob = 1) ----
    Factor("B_v",      (BAS,), "v",     0.432, 0.72, source=V_SOURCE),
    Factor("B_b",      (BAS,), "b",     0.55, 1.35,
           source="Verzweigungsrate / Spitzendichte (Oyarte Galvez et al. 2025, Fig. b): ≈ 0,044 h⁻¹ ≈ 1,05 /d am Maximum, "
                  "≈ 0,65 /d an der Flanke; Band 2 SEM ≈ ±40 %. Zählt dort alle Verzweigungen (auch Runner-Laterale)"),
    Factor("B_hlt",    (BAS,), "hlt",   5.0, 9.0,
           source="BAS degenerieren nach ~7 d (Bago 1998; Kokkoris 2026)"),
    Factor("B_theta",  (BAS,), "theta", math.radians(50), math.radians(80),
           source="wie R_theta; Messung in Oyarte Galvez 2025 mit 881 BAS-Verzweigungspunkten"),

    # ---- gemeinsam für Runner und BAS ----
    Factor("distTH", (RUNNER, BAS), "distTH", *_pm(PH, 0.01)),   # nominal 0,01 cm = dx
    Factor("ana",    (RUNNER, BAS), "ana",    0.5, 1.0),        # Wahrscheinlichkeit, nominal 1 = Rand
    Factor("v_cv",   (RUNNER, BAS), "v_cv",   0.15, 0.40,
           source="Streuung der gemessenen Spitzengeschwindigkeiten (Oyarte Galvez et al. 2025), nominal CV ≈ 0,29"),
]

# Parameter, die NIE variiert werden (nur zur Dokumentation / Prüfung)
FIXED = ["successor percentage (0.75 BAS, 0.25 Runner)", "la", "lb", "BAS lmax (≈ 3 × dx)", "a",
         "tropismT", "tropismN", "tropismS", "dx", "dxMin", "ldelay", "b_prob"]

# --------------------------------------------------------------------------
# Simulation (wie petridishsimauto.py)
# --------------------------------------------------------------------------
N_RINGS = 15
HYPHAE_LAYER = 0.1          # Hyphen wachsen auf dem Cellophan, dünne Schicht oben [cm]
ROOT_DAYS = 5               # Wurzelphase [d]
STEPS_PER_DAY = 24          # dt = 1 h
HOURS_HYPHAE = 60           # Simulation nach t* [h], wie die Darstellungen im Paper
MAX_WAIT_DAYS = 30          # ohne t* nach so vielen Tagen: Lauf zählt als fehlgeschlagen
# Kein Abbruch bei stagnierender Front: in der SA ist das ein Ergebnis (kleines v), kein Fehler
MIN_GROWING_TIPS = None
STALL_FRACTION = None
SAVE_VTP = False            # finale Geometrie als .vtp speichern (für Netzwerkauswertungen, groß)

# --------------------------------------------------------------------------
# Outputs (evaluate.py)
# --------------------------------------------------------------------------
SNAPSHOT_HOURS = [10, 30, 60]   # Auswertung von Mengen nach t* [h], wie im Paper (0, 10, 30, 60 h)
# Mengen zu jedem Zeitpunkt in SNAPSHOT_HOURS
SNAPSHOT_OUTPUTS = ["total_length", "active_tips", "bas_share", "anastomoses", "front_position"]
# Wellengrößen über das ganze Fenster bis max(SNAPSHOT_HOURS), stehen beim letzten Zeitpunkt
WAVE_OUTPUTS = ["front_velocity", "K1", "K2", "wave_width", "encounter_C"]
# Gütekriterien der Wellen-Fits wie petridishstats.py
R2_MIN_LENGTH = 0.9
R2_MIN_TIPS = 0.5
SATURATION_MIN = 0.8
MIN_RINGS = 3
ENCOUNTER_MIN_TIPS = 5

# --------------------------------------------------------------------------
# Phase 0
# --------------------------------------------------------------------------
PHASE0_NOISE_SEEDS = 50
PHASE0_DX_VALUES = [0.01, 0.005, 0.0025]   # nominal, halb, viertel (Laufzeit steigt etwa mit 1/dx!)
PHASE0_DX_SEEDS = 10
PHASE0_RANGE_POINTS = 20    # zufällige Punkte im Faktorraum, je ein Lauf: Abbruchquote und Laufzeit am Rand

# --------------------------------------------------------------------------
# Phase 1 (Morris)
# --------------------------------------------------------------------------
MORRIS_LEVELS = 4            # p
MORRIS_TRAJECTORIES = 20     # r (aus MORRIS_CANDIDATES nach Spreizung ausgewählt)
MORRIS_CANDIDATES = 300
MORRIS_REPLICATES = 5        # Phase 0: Seed-CV 0,3–0,97 je Output -> 5 statt 3
BASE_SEED = 20261021
