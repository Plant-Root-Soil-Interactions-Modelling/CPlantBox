# Sensitivitätsanalyse Hyphenmodell — Phase 0 + Morris

Gerüst für Rauschtest, `dx`-Konvergenz und Morris-Screening auf SLURM.
Mit einem Spielzeugmodell (`dummy_model.py`) komplett durchgetestet; für die echte
Simulation muss nur `model_hook.py` angeschlossen werden.

## Was du anpassen musst

| Datei | Was |
|---|---|
| `template/hyphae_template.xml` | durch dein vollständiges Nominal-Parameterfile ersetzen |
| `model_hook.py` | `simulate()` an eure Python-Bindings anschließen (Seed setzen, Snapshots speichern), `carbon_cost()` mit eurer Definition |
| `config.py` | **Faktorbereiche** (aktuell alle ±25 % Platzhalter), `SIM_TIME`, `SNAPSHOT_TIMES`, `DOMAIN`, `REGIONS` |
| `slurm/*.sbatch` | Umgebung laden (module/venv), ggf. Partition/Account, Speicher |

Falls eure bestehende Auswertung die Outputs schon berechnet: in `evaluate.py`
`region_outputs()` durch Aufrufe eurer Funktionen ersetzen. Wichtig ist nur das
Ausgabeformat (eine Zeile pro `time, region, output, value`).

## Ablauf

```bash
python3 check_config.py                      # Bereiche, Nominalwerte, Laufanzahl, CPU-Stunden

# ---- Phase 0: Rauschen + dx-Konvergenz (80 Läufe) ----
python3 design.py phase0
bash slurm/submit.sh results/design_phase0.csv
python3 analyze_phase0.py                    # nach Abschluss

# ---- Phase 1: Morris (r × (k+1) × R Läufe) ----
python3 design.py morris
bash slurm/submit.sh results/design_morris.csv
python3 analyze_morris.py
```

Lokal testen ohne Cluster und ohne echtes Modell:

```bash
export SA_DUMMY=1
python3 run_one.py results/design_phase0.csv 0-4
python3 evaluate.py results/design_phase0.csv 0-4
```

**Abgebrochene Läufe** (Timeout, Knotenausfall): `submit.sh` einfach erneut aufrufen,
fertige Läufe werden übersprungen. `results/status_*.csv` zeigt, was fehlt.

**Neue Outputs / andere Teilbereiche / andere Zeitpunkte:** `evaluate.py` bzw.
`config.REGIONS` ändern, dann `bash slurm/submit.sh <design> --eval-only`.
Keine Neusimulation nötig, solange die Snapshots die Rohdaten enthalten.

## Designentscheidungen

- **Variiert werden die `value`s**, nicht die `dev`s. Mit `DEV_MODE = "relative"`
  wird `dev` proportional mitskaliert (konstanter CV), sonst entstehen bei kleinen
  Mittelwerten massenhaft negative Ziehungen.
- **b_prob als Struktur, nicht als Faktor:** Puller ohne `b`, Runner ohne `ln/la/lb`.
- **Kontrollfaktor `P_b_ctrl`** (b der Puller, wirkungslos bei b_prob = 0): muss μ* = 0
  liefern. Sonst greift der Schalter nicht so wie gedacht.
- **Successor-Aufteilung** als ein Faktor p, geschrieben als `"p, 1-p"`.
- **`dx` nur in Phase 0** (Konvergenztest), nie in der Sensitivitätsanalyse.
- **Gemeinsame Seeds entlang einer Trajektorie** (Common Random Numbers): Elementary
  Effects sind Differenzen, gleiche Seeds nehmen einen Großteil des Rauschens raus.
- **Rauschniveau aus den Replikaten:** Jedes Replikat ist eine eigene konsistente
  Trajektorie. Die Streuung der EEs zwischen Replikaten ergibt pro Faktor das μ*,
  das reines Seed-Rauschen erzeugen würde. Faktoren mit μ* ≤ 2 × diesem Niveau gelten als
  „im Rauschen" (Spalte `in_noise`, schwarze Markierung im Ranking-Plot).
- **μ* in normierten Einheiten:** Output-Änderung pro Durchlaufen des Faktorbereichs,
  zwischen Faktoren vergleichbar. `mu_star_rel` normiert zusätzlich je Output.
- Morris ist selbst implementiert (`morris.py`, Trajektorienauswahl nach Campolongo
  et al. 2007), damit Replikat- und Seed-Logik sauber passen. Ergebnisse entsprechen
  `SALib.analyze.morris`.

## Was Phase 0 entscheidet

- `phase0_noise.csv`: CV je Output. Bei hohem CV (> ~0.3) eher `MORRIS_REPLICATES` erhöhen.
  Spalte `frac_zero` zeigt spärlich besetzte Teilbereiche, dort sind varianzbasierte
  Maße später mit Vorsicht zu lesen.
- `phase0_dx.csv`: Weicht nominales `dx` signifikant (> 2 %) vom feinsten ab, ist
  `dx` zu grob. Erst klären, dann Morris starten.
- Laufzeit und Speicher aus `status_phase0.csv` bzw. SLURM-Logs: `--time`/`--mem` anpassen.

## Ergebnisse

- `results/morris_indices.csv`: μ*, KI, μ, σ, Rauschniveau je `time × region × output × factor`
- `results/plots/morris_ranking_t*.png`: μ* mit Bootstrap-KI, gesamter Zielbereich
- `results/plots/morris_musigma_t*.png`: μ* vs. σ (σ ≳ μ*: nichtlinear/Interaktionen)
- `results/plots/morris_regions_t*.png`: relative Sensitivität je Teilbereich

`beispiel_dummy/` enthält diese Plots aus dem Testlauf mit dem Spielzeugmodell.
Sie zeigen nur das Format, inhaltlich bedeuten sie nichts.

## Literatur

- Morris (1991), *Technometrics* 33:161–174
- Campolongo, Cariboni, Saltelli (2007), *Environ. Model. Softw.* 22:1509–1518 (μ*, Trajektorienauswahl)
