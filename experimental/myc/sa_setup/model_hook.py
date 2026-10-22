"""
SCHNITTSTELLE ZU DEINER SIMULATION — die einzige Datei, die du anpassen musst.

Zwei Funktionen:

1. simulate(xml_path, seed, outdir) -> list[Path]
   Startet die Simulation über eure Python-Bindings mit dem gegebenen Parameterfile
   und Seed und speichert zu jedem Zeitpunkt in config.SNAPSHOT_TIMES die ROHAUSGABE
   (so vollständig wie möglich, damit spätere Auswertungen ohne Neusimulation gehen).

2. load_snapshot(path) -> dict
   Liest eine gespeicherte Rohausgabe und gibt sie in einem einheitlichen Format zurück:
       "seg_start": (N, 3)  Startpunkte der Segmente [cm]
       "seg_end":   (N, 3)  Endpunkte [cm]
       "radius":    (N,)    Radius [cm]
       "subtype":   (N,)    subType der zugehörigen Hyphe
       "tips":      (M, 3)  Positionen der aktiven Spitzen
       "anastomoses": (A, 3) Positionen der Anastomosen
   Wenn eure bestehende Auswertung andere Größen braucht, einfach Einträge ergänzen
   und in evaluate.py verwenden.

3. carbon_cost(snap) -> float  (pro Region aufgerufen, s. evaluate.py)

Solange die Umgebungsvariable SA_DUMMY=1 gesetzt ist, wird ein Spielzeugmodell
verwendet (dummy_model.py), damit die Pipeline ohne die echte Simulation testbar ist.
"""
import os
from pathlib import Path

import numpy as np

import config

USE_DUMMY = os.environ.get("SA_DUMMY") == "1"


def simulate(xml_path: Path, seed: int, outdir: Path) -> list[Path]:
    if USE_DUMMY:
        import dummy_model
        return dummy_model.simulate(xml_path, seed, outdir)

    # ------------------------------------------------------------------
    # TODO: an eure Bindings anschließen. Skizze (CPlantBox-artig):
    #
    #   import plantbox as pb
    #   plant = pb.MycorrhizalPlant()          # bzw. eure Klasse
    #   plant.readParameters(str(xml_path))
    #   plant.setSeed(seed)                    # WICHTIG: Seed wirklich setzen
    #   plant.initialize()
    #   paths, t = [], 0.0
    #   for t_snap in config.SNAPSHOT_TIMES:
    #       plant.simulate(t_snap - t)
    #       t = t_snap
    #       p = outdir / f"snapshot_t{t_snap:g}.npz"
    #       np.savez_compressed(p, seg_start=..., seg_end=..., radius=...,
    #                           subtype=..., tips=..., anastomoses=...)
    #       paths.append(p)
    #   return paths
    # ------------------------------------------------------------------
    raise NotImplementedError("model_hook.simulate an die Simulation anschließen")


def load_snapshot(path: Path) -> dict:
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


# Kohlenstoffkosten — PLATZHALTER: Volumen × Konstante. Durch eure Definition ersetzen.
CARBON_PER_CM3 = 1.0   # TODO: echte Umrechnung (z. B. g C pro cm³ Hyphenvolumen)


def carbon_cost(seg_len, radius, subtype) -> float:
    return float(CARBON_PER_CM3 * np.sum(np.pi * radius**2 * seg_len))
