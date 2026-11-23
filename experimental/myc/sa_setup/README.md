# Sensitivitätsanalyse Hyphenmodell — Phase 0 + Morris

Gerüst für Rauschtest, `dx`-Konvergenz und Morris-Screening auf SLURM, für das
Petrischalen-Setup von `petridishsimauto.py` (Oyarte Galvez et al. 2025) mit
`tomatoparameters/TwoHyphaePlusBAS.xml` als Nominalmodell. Mit einem Spielzeugmodell
(`dummy_model.py`) komplett durchgetestet.

## Was du anpassen musst

| Datei | Was |
|---|---|
| `config.py` | **Faktorbereiche** (`distTH` noch ±25 % Platzhalter), ggf. `HOURS_HYPHAE`, `SNAPSHOT_HOURS` |
| `slurm/*.sbatch` | aktivieren `~/cpbenv` (braucht zusätzlich SALib, scikit-learn, networkx), Partition `cpu256`; ggf. Zeit/Speicher |

`model_hook.py` ruft `AMFSimulationSetup.makeSimulation` mit dem Parameterfile des Laufs auf
und speichert die Größen je Ring und Stunde wie `petridishsimauto.py` (`raw.npz`).
Änderungen am Nominalmodell gehören in `TwoHyphaePlusBAS.xml`, nicht hierher.

## Ablauf

Die Jobs werden wie die Petrischalen-Läufe aus dem **Home-Verzeichnis** eingereicht, und
**jedes Python-Skript läuft als Job** (nicht auf dem Login-Knoten). Einmalig:

```bash
cd ~
cp CPlantBox/experimental/myc/sa_setup/slurm/*.sbatch CPlantBox/experimental/myc/sa_setup/slurm/submit_smoke.sh .
cp CPlantBox/experimental/myc/sa_setup/slurm/submit.sh submit_sa.sh
source cpbenv/bin/activate && pip install SALib scikit-learn networkx
mkdir -p slurmOut slurmErr
```

Alles aus `~`. Einzelne Python-Schritte mit `sbatch sa_run.sbatch <skript> [argumente]`,
ihre Ausgabe steht in `slurmOut/sa_run-<job>.out`:

```bash
cd ~
sbatch sa_run.sbatch check_config.py         # Bereiche, Nominalwerte, Laufanzahl, CPU-Stunden

# ---- Smoke-Test: 3 echte Läufe (nominal, alle Faktoren unten, alle oben), danach Bericht ----
bash submit_smoke.sh                         # oder: bash submit_smoke.sh --quick (12 h nach t*)
# Bericht: slurmOut/sa_run-<job>.out -> alles OK? Laufzeit und Speicher -> sim.sbatch

# ---- Phase 0: Rauschen, dx-Konvergenz, Stichprobe des Faktorraums (100 Läufe) ----
sbatch sa_run.sbatch design.py phase0        # warten, bis der Job fertig ist
bash submit_sa.sh results/design_phase0.csv  # Designpfad relativ zu sa_setup
# Sammeln und analyze_phase0.py laufen am Ende automatisch: slurmOut/sa_collect-<job>.out

# ---- Test: Nominalpunkt über 120 h (10 Läufe), sind die Wellen-Fits dann oft genug definiert? ----
sbatch sa_run.sbatch design.py long          # warten, bis der Job fertig ist
SIMTIME=12:00:00 bash submit_sa.sh results/design_long.csv
# Ergebnis: slurmOut/sa_collect-<job>.out (analyze_long.py), Fenster 60, 90 und 120 h

# ---- Phase 1: Morris (r × (k+1) × R Läufe) ----
sbatch sa_run.sbatch design.py morris        # warten, bis der Job fertig ist
bash submit_sa.sh results/design_morris.csv
# analyze_morris.py läuft am Ende automatisch, Plots in sa_setup/results/plots/
```

Ein Schritt noch einmal von Hand, z. B. nach Änderungen an der Analyse:
`sbatch sa_run.sbatch analyze_morris.py`.

Liegt CPlantBox woanders, `SA_DIR=/pfad/zu/sa_setup` vor `sbatch` bzw. `bash submit_sa.sh` setzen.
Nach Änderungen an den sbatch-Dateien im Repo die Kopien im Home-Verzeichnis erneuern.

Lokal testen ohne Cluster und ohne echtes Modell:

```bash
export SA_DUMMY=1
python3 design.py phase0
python3 run_one.py results/design_phase0.csv 0-4
python3 evaluate.py results/design_phase0.csv 0-4
```

**Abgebrochene Läufe** (Timeout, Knotenausfall): `submit_sa.sh` einfach erneut aufrufen,
fertige Läufe werden übersprungen. `results/status_*.csv` zeigt, was fehlt. Läufe ohne t*
nach `MAX_WAIT_DAYS` zählen als fehlgeschlagen (Fehlermeldung in `runs/<id>/meta.json`).

**Neue Outputs / andere Zeitpunkte:** `evaluate.py` bzw. `config.SNAPSHOT_HOURS` ändern, dann
`cd ~ && bash submit_sa.sh <design> --eval-only`. Keine Neusimulation nötig, `raw.npz` enthält
alle Größen je Ring und Stunde. Für Netzwerkgrößen (`AMFNetwork.py`) braucht es die
Geometrie: dafür `SAVE_VTP = True` setzen.

## Faktoren (13)

| Gruppe | Faktoren | Bereich |
|---|---|---|
| Runner (subType 1) | `R_v`, `R_hlt`, `R_lmax`, `R_ln`, `R_theta` | v 0,432–0,72 cm/d (180–300 µm/h), hlt 4–10 d, lmax 1,5–5 cm, ln 0,1–0,6 cm (log), θ 50–80° |
| BAS (subType 2) | `B_v`, `B_b`, `B_hlt`, `B_theta` | v wie Runner, b 0,55–1,35 /d, hlt 5–9 d (nominal 7 d), θ 50–80° |
| gemeinsam | `distTH`, `ana`, `v_cv` | ±25 %, 0,5–1, CV 0,15–0,40 |

Die Quellen der Bereiche stehen je Faktor in `config.py` (`source`).
| Kontrolle | `R_b_ctrl` | b der Runner, wirkungslos bei b_prob = 0 |

Fest: Successor-Aufteilung (0,75 BAS / 0,25 Runner), la, lb, BAS-lmax (≈ 3 × dx),
Radius a, Tropismus, dx, b_prob.

## Outputs

| Art | Outputs | Zeit |
|---|---|---|
| Mengen | `total_length` [cm], `active_tips`, `bas_share`, `anastomoses` (seit t*), `front_position` (spitzengewichteter Radius, mm) | 10, 30, 60 h nach t* |
| Welle | `front_velocity` [mm/d], `K1` [µm mm⁻²], `K2` [mm⁻²], `wave_width` [h], `encounter_C`, `beta` [µm/h] (a = β ρ n, Ziel ≈ 23) | Fits über 0–60 h, stehen bei 60 h |

Die Wellengrößen kommen aus den Fits von `AMFDensities` mit den Gütekriterien von
`petridishstats.py` (R², Sättigung, mind. 3 Ringe). Scheitern die Fits, fehlt der Wert,
und der Lauf zählt für diesen Output wie ein fehlgeschlagener. Die Anzahl Ringe hinter
den Fits steht als Diagnose (`n_wave_rings`, `n_K1_rings`, `n_tip_rings`) in `outputs.csv`,
ebenso `wait_days` (Zeit bis t*). Diagnosen gehen nicht in die SA ein.

## Designentscheidungen

- **Variiert werden die `value`s.** Mit `DEV_MODE = "relative"` wird `dev` proportional
  mitskaliert (konstanter CV). Ausnahme `v_cv`: setzt dev von v auf CV × v, für beide Subtypen.
- **Gemeinsame Faktoren** (`distTH`, `ana`, `v_cv`) werden auf Runner und BAS geschrieben.
- **b_prob als Struktur, nicht als Faktor.** b wird nur beim Tip-Splitting gelesen
  (`Hyphae.cpp`), daher ist `R_b_ctrl` ein **Kontrollfaktor**: muss μ* = 0 liefern.
- **Kein Abbruch bei stagnierender Front** (`MIN_GROWING_TIPS = None`): in der SA ist eine
  langsame Front ein Ergebnis, kein Fehler. Nur ohne t* wird abgebrochen.
- **`dx` nur in Phase 0** (Konvergenztest), nie in der Sensitivitätsanalyse.
- **Gemeinsame Seeds entlang einer Trajektorie** (Common Random Numbers): Elementary
  Effects sind Differenzen, gleiche Seeds nehmen einen Großteil des Rauschens raus.
- **Rauschniveau aus den Replikaten:** Jedes Replikat ist eine eigene konsistente
  Trajektorie. Die Streuung der EEs zwischen Replikaten ergibt pro Faktor das μ*,
  das reines Seed-Rauschen erzeugen würde. Faktoren mit μ* ≤ 2 × diesem Niveau gelten als
  „im Rauschen" (Spalte `in_noise`, schwarze Markierung im Ranking-Plot).
- **μ* in normierten Einheiten:** Output-Änderung pro Durchlaufen des Faktorbereichs,
  zwischen Faktoren vergleichbar. `mu_star_rel` normiert zusätzlich je Output.
- **Morris über SALib** (`pip install SALib`): Design mit `SALib.sample.morris`
  (Trajektorienauswahl nach maximaler Spreizung, `local_optimization`), Auswertung mit
  `SALib.analyze.morris` auf dem Replikat-Mittel. Faktoren werden auf `[0,1]` gesampelt
  und in `design.to_physical` umgerechnet (auch log-uniform). Selbst berechnet wird nur
  das Rauschniveau, weil SALib die einzelnen EEs nicht herausgibt.
- **Fehlgeschlagene Läufe:** Ein Replikat einer Trajektorie zählt nur, wenn alle k+1
  Punkte vorliegen; eine Trajektorie zählt, wenn mind. ein Replikat vollständig ist.

## Was Phase 0 entscheidet

- `phase0_noise.csv`: CV je Output. Bei hohem CV (> ~0.3) eher `MORRIS_REPLICATES` erhöhen.
  Spalte `frac_zero` zeigt spärlich besetzte Outputs.
- `phase0_dx.csv`: Weicht nominales `dx` signifikant (> 2 %) vom feinsten ab, ist
  `dx` zu grob. Besonders für Anastomosen (distTH = dx) und BAS (lmax ≈ 3 × dx) prüfen.
  Kleines dx kostet viel Laufzeit: `--time` in `sim.sbatch` dafür erhöhen.
- `phase0_status.csv`: Wie viele der 20 zufälligen Punkte im Faktorraum t* erreichen, und wie
  lange sie laufen. Viele Abbrüche: Bereiche zu breit, sonst fehlen bei Morris ganze Trajektorien.
- Laufzeit und Speicher: `--time`/`--mem` in den sbatch-Dateien anpassen.

## Ergebnisse

- `results/morris_indices.csv`: μ*, KI, μ, σ, Rauschniveau je `time × region × output × factor`
- `results/plots/morris_ranking_t*.png`: μ* mit Bootstrap-KI (SALib: μ* ± Halbbreite) je Output
- `results/plots/morris_musigma_t*.png`: μ* vs. σ (σ ≳ μ*: nichtlinear/Interaktionen)
- `results/plots/morris_overview.png`: relative Sensitivität aller Outputs und Zeitpunkte

`beispiel_dummy/` enthält Plots aus einem Testlauf mit dem Spielzeugmodell.
Sie zeigen nur das Format, inhaltlich bedeuten sie nichts.

## Kalibrierung (`calibrate.py`)

Nullstellensuche `F(x) = 0` mit `F_i = log(y_i / y*_i) / log(1 + tol_i)` über die Zielwerte
β ≈ 23 µm/h, K1 ≈ 1400 µm mm⁻², BAS-Anteil 0,3 und K2 ≈ 0,15 mm⁻². Ein Ziel gilt als
erreicht, wenn |F_i| ≤ 1 ist. Das Verfahren ist Broyden (mehrdimensionale Sekante) im normierten
Parameterraum. Die Start-Jacobimatrix kommt aus großen Vorwärtsdifferenzen. Alle
Auswertungen verwenden dieselben Seeds (CRN), die Outputs sind über die Seeds gepoolt, der
Standardfehler kommt aus einem Bootstrap. Parameter, Suchbereiche, Ziele und Toleranzen
stehen oben in `calibrate.py` (`PARAMS`, `TARGETS`).

```bash
cd ~ && cp CPlantBox/experimental/myc/sa_setup/slurm/calibrate.sbatch .
sbatch sa_run.sbatch calibrate.py --check    # Bereiche, Startpunkt, Anzahl Auswertungen
sbatch calibrate.sbatch                      # 16 Kerne, 16 Seeds je Auswertung
```

Ergebnisse in `calibration/`: `evals.csv` (jede Auswertung mit Werten, SE, F),
`jacobian.csv` (geschätzte Sensitivitäten dF/du, zeigt, welcher Parameter welches Ziel treibt),
`best_params.xml`. Neustartfähig: denselben Befehl erneut einreichen.
Lokal testen: `SA_DUMMY=1 python3 calibrate.py --seeds 8`.

## Literatur

- Morris (1991), *Technometrics* 33:161–174
- Campolongo, Cariboni, Saltelli (2007), *Environ. Model. Softw.* 22:1509–1518 (μ*, Trajektorienauswahl)
- Oyarte Galvez et al. (2025): Versuchsaufbau, Spitzengeschwindigkeit 240 µm/h
