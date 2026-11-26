# Anleitung: Läufe auf dem Cluster

#amf #cluster #sensitivitaetsanalyse #kalibrierung

Reihenfolge nach B14: **Phase 0 + Langzeittest + Inokulum-Morris → auswerten → Morris → auswerten, Parameter wählen →
Kalibrierung → Validierung**.

**Automatischer Ablauf (Standard):** `submit_all.sh` reicht Smoke-Test, Phase 0 und Langzeittest ein. Danach prüft ein
Gate (`phase0_gate.py`) die Ergebnisse automatisch:
- **GO:** Morris und die Kalibrierung starten von selbst. Die Kalibrierung läuft mit allen Parametern mit `calibrate=True`
  und mit `between_cv` aus Phase 0.
- **STOP:** Es wird nichts weiter gestartet.

In beiden Fällen kommt eine E-Mail mit dem Bericht (siehe Abschnitt 1a). Die Auswahl der Parameter nach Morris (Abschnitt 4)
bleibt Handarbeit; mit ihr kann die Kalibrierung danach noch einmal laufen.

Alle Befehle aus dem **Home-Verzeichnis** (`cd ~`). Ausgaben der Jobs landen in `~/slurmOut`, Fehler in `~/slurmErr`,
Ergebnisse in `~/CPlantBox/experimental/myc/sa_setup/results/` bzw. `.../calibration/`. Python nie auf dem Login-Knoten
starten, sondern über `sa_run.sbatch`.

Abkürzung für diese Anleitung:
```bash
SA=~/CPlantBox/experimental/myc/sa_setup
S=$SA/slurm
```

---

## 0. Vorbereitung (einmalig)

**Code auf den Cluster bringen** (lokal pushen, auf dem Cluster pullen):
```bash
cd ~/CPlantBox && git fetch && git checkout AMF_Anastomosis && git pull
git log -1 --oneline      # muss den Commit "sa_setup: rework of sensitivity analysis and calibration" (oder neuer) zeigen
```
Am C++-Code wurde nichts geändert, `plantbox` muss nicht neu gebaut werden.

**Python-Umgebung prüfen** (SALib, scikit-learn, networkx, scipy, pandas, matplotlib):
```bash
source ~/cpbenv/bin/activate && pip install SALib scikit-learn networkx
mkdir -p ~/slurmOut ~/slurmErr
```

**Alte Ergebnisse sichern.** Das Nominalmodell hat sich geändert; neue Designs überschreiben `results/design_*.csv`,
`outputs_*.csv` und die Kalibrierung schreibt `calibration/evals.csv` neu. Die alten Läufe in `runs/` bleiben (andere IDs).
```bash
cd $SA
[ -d results ] && mv results results_v1
[ -d calibration ] && mv calibration calibration_v1
cd ~
```

**Konfiguration prüfen:**
```bash
sbatch $S/sa_run.sbatch check_config.py          # Faktoren, Bereiche, Priors, Zahl der Läufe, CPU-Stunden
sbatch $S/sa_run.sbatch calibrate.py --check     # 19 Parameter, 21 Ziele mit Quellen, Kosten pro Iteration
```
Ausgabe in `~/slurmOut/sa_run-<job>.out`. Keine Zeile mit `ERROR`.

---

## 1a. Alles automatisch mit Gate und E-Mail

```bash
cd ~ && bash $S/submit_all.sh --mail <adresse>        # oder: export SA_MAIL=<adresse>
```

Ablauf:
1. Smoke-Test und Smoke-Gate.
2. Parallel: Phase 0 (210 Läufe) und Langzeittest (20 Läufe, 150 h).
3. Sind beide gesammelt und ausgewertet, läuft das **Phase-0-Gate** (`slurm/gate.sbatch` → `phase0_gate.py`):
   - **GO:** Es startet die Morris-Stufe und die Kalibrierung (`--cal-args "..."` wird weitergegeben).
   - **STOP:** Nichts weiter.
   - Bericht in `results/gate_report.txt`, E-Mail mit dem Bericht als Text.
   - Fallback ohne `mail`/`sendmail` auf dem Knoten: SLURM-Mail mit dem Ergebnis im Betreff; der Bericht liegt dann nur in
     der Datei.
   - Bricht eine Stufe ab (z. B. das Smoke-Gate scheitert), kommt ebenfalls eine STOP-Mail.

Optionen:
- `--no-long`: Das Gate prüft nur Phase 0.
- `--no-morris` / `--no-calibration`: Das Gate startet das jeweilige nicht.

**STOP** (die folgenden Stufen würden kaputte oder unbrauchbare Läufe rechnen). Schwellen in `phase0_gate.py` (`GATE_*`):

| Prüfung | STOP, wenn |
|---|---|
| Rausch-Läufe (Nominalpunkt) mit Ergebnis | < 80 % |
| Zufallspunkte im Faktorraum mit Ergebnis | < 50 % (WARN < 80 %) |
| Laufzeit der Zufallspunkte, 90 %-Quantil | > 75 % des Zeitlimits von `sim.sbatch` (8 h) |
| Wellenziele (alpha, beta, K1, K2, front_velocity, front_tip_speed) pro Rausch-Lauf definiert | < 30 % (WARN < 60 %) |
| dx 0.01 gegen 0.005 cm, jedes Kalibrierziel | \|Abweichung\| > 10 % und CI ohne 0 (WARN: CI außerhalb ±2 %) |
| Ringdiagnostik | fehlt ganz (andere Ring-Hinweise: WARN) |
| Langzeittest (Wurzeln aus) mit Ergebnis | < 70 % |
| Wellenziele 120 → 150 h, Median der gepaarten Änderung | > ±25 % (WARN > ±10 %): Fenster `HOURS_HYPHAE` zu kurz |
| Datei fehlt (Sammeln oder Analyse gescheitert) | immer |

Nur Information im Bericht: dt, Schichtdicke, CRN, Laufzeiten, `between_cv`.

`results/between_cv.json`:
- Robuste Streuung zwischen den Seeds (1.4826 · MAD / |Median|) der Ziele aus einer Kolonie.
- `calibrate.py` liest die Datei bei jedem Start. Ein von Hand in `calibrate.TARGETS` gesetzter Wert (≠ 0) hat Vorrang.

**Nach einem STOP:**
1. Ursache beheben (Konfiguration, Zeitlimit).
2. Neu einreichen mit `bash $S/submit_all.sh --no-smoke --mail <adresse>`. Fertige Läufe werden übersprungen, also wird
   nur Fehlendes gerechnet.

Soll es trotz WARN oder STOP weitergehen, startet man von Hand: Morris mit Abschnitt 3, die Kalibrierung mit Abschnitt 5.

---

## 1. Stufe 1 von Hand: Smoke-Test, Phase 0, Langzeittest, Inokulum-Morris

```bash
cd ~ && bash $S/submit_all.sh --no-morris --no-calibration
```
Das reicht ein (mit SLURM-Abhängigkeiten):
1. **Smoke-Test** (3 Läufe: nominal, alle Faktoren unten, alle oben) und das **Gate**: Läuft der Nominalpunkt nicht sauber
   durch, bricht SLURM alles Folgende ab. Ausgabe: `slurmOut/sa_run-<gate>.out` (OK / PROBLEMS pro Punkt, Laufzeit, Speicher).
2. danach parallel:
   - **Phase 0**: 210 Läufe (Rauschen 50, dx 20, Zufallspunkte 20, dt 30, Schichtdicke 30, anaMinDist 40, CRN 20)
   - **Langzeittest**: 20 Läufe über 150 h (Wurzeln aus / an)
   - **Inokulum-Morris**: 180 Läufe

**Zeitlimits:**
- Simuliert werden 120 h nach t* (`HOURS_HYPHAE`), ein Nominallauf dauert ~1.6 h. `sim.sbatch` hat 8 h pro Lauf.
  Die dt-Läufe mit 0.25 h brauchen etwa die 4-fache Zeit eines Nominallaufs, also eventuell mehr: Zeigt der Smoke-Test
  mehr als ~2 h für den Nominallauf, Phase 0 mit mehr Zeit nachreichen (siehe Abschnitt 7).
- Langzeittest (150 h, nominal ~2.5 h pro Lauf; die Laufzeit wächst etwa quadratisch mit der simulierten Zeit): Standard
  16 h (`SIMTIME_LONG`). Mit weiterwachsenden Wurzeln entstehen mehr Hyphen; bei Zeitüberschreitungen großzügiger:
  ```bash
  SIMTIME_LONG=24:00:00 bash $S/submit_all.sh --no-morris --no-calibration
  ```

**Beobachten:**
```bash
squeue -u $USER                                     # laufende / wartende Jobs
cat $SA/results/status_phase0.csv | cut -d, -f2 | sort | uniq -c    # ok / excluded / failed / missing
grep -h -i -E 'TIME LIMIT|oom' ~/slurmErr/*.err | sort | uniq -c   # Abbrüche wegen Zeit oder Speicher
```

---

## 2. Stufe 1 auswerten

Die Analysen laufen automatisch nach dem Sammeln (`slurmOut/sa_collect-<job>.out`). Von Hand neu:
```bash
sbatch $S/sa_run.sbatch analyze_phase0.py
sbatch $S/sa_run.sbatch analyze_long.py
sbatch $S/sa_run.sbatch analyze_morris.py inoculum
sbatch $S/sa_run.sbatch analyze_rings.py            # Ringauflösung 15 vs. 5 auf den Rausch-Läufen
```

| Datei (in `results/`) | Frage | Folge |
|---|---|---|
| `phase0_warnings.txt` | „ATTENTION …“? Wenige Kreuzungen bis t*, Spitzen jenseits 32 mm | Eintrittspunkt länger mitteln bzw. `RING_RADIUS` vergrößern |
| `phase0_noise.csv` | CV pro Output | > ~0.3: mehr Replikate (`MORRIS_REPLICATES`) |
| Ausgabe „Spread between seeds …“, `between_cv.json` | CV der Ziele aus einer Kolonie | schreibt das Gate automatisch (`between_cv.json`, robust); nur bei Bedarf von Hand in `calibrate.TARGETS` überschreiben |
| `phase0_dx.csv` | Abweichung 0.01 gegen 0.005 cm > 2 % (und größer als das Seed-Rauschen)? | gezielter Zusatztest bei 0.0025 cm mit wenigen Seeds, ggf. dx verkleinern |
| `phase0_dt.csv` | Rest-Abhängigkeit von dt (dt = 1 h ist festgelegt) | nur dokumentieren; Anastomosen pro Spitze sollten gleich sein |
| `phase0_layer.csv` | Wie stark reagieren alpha, Anastomosen auf die Schichtdicke (0.02 nominal vs. 0.005 / 0.1 cm)? | Vorschau auf den Morris-Effekt von `hyphae_layer`; Nominalwert halten oder `HYPHAE_LAYER` (und den Nominalwert des Faktors) ändern |
| `phase0_anaL.csv` | Wirkung der Mindestpfadlänge (×1 / 1.5 / 6 gegen nominal ×3) | `ana_L_scale` in der Kalibrierung lassen oder fixieren |
| `phase0_crn.csv` | Verhältnis gleiche/unabhängige Seeds deutlich < 1? | sonst mehr Seeds in Kalibrierung und Morris |
| `phase0_bas_share.csv` | BAS-Anteil Subtyp vs. Bildkriterien | Einordnung C15 |
| `phase0_status.csv` | Abbrüche an den Rändern des Faktorraums (Gruppe range) | viele: Bereiche prüfen, `MAX_WAIT_DAYS` |
| `phase0_noise.csv`: `root_efficiency_rel`, `global_efficiency_rel`, `cost_rel` | zwischen 0 und 1, Ĉ > 0? | sonst Netz nicht zusammenhängend |
| `phase0_noise.csv`: `branch_angle_*`, `persistence_length`, `front_tip_speed`, `wave_width` | Abstand zum Experiment, Rauschen | `SPEED_FROM`; `wave_width` als Ziel aufnehmen? |
| `long_roots_compare.csv`, `plots/long_influx.png` | Unterschied Wurzeln an/aus innerhalb 120 h? | ja: Wurzeln weiterwachsen lassen (`turn_off_root_growth`) |
| `long_summary.csv` | K1/K2 bei 120 h definiert und gegenüber 150 h konvergiert? | sonst `HOURS_HYPHAE` erhöhen |
| `inoculum_vs_hyphae.csv` | Verhältnis ≪ 1? | Normierung auf t* nimmt die Wurzelphase heraus |
| `rings_compare_phase0.csv` | Hängen K1, K2, Wellenfits an der Ringbreite? | Einordnung |

**Änderungen aus der Auswertung** (Konfiguration, `between_cv`, Numerik) lokal machen, committen, auf dem Cluster pullen.
Ändert sich dabei das Modell (dx, dt, Schichtdicke, Ringe), muss Phase 0 nicht komplett neu laufen, aber die folgenden
Stufen mit dem neuen Stand.

---

## 3. Stufe 2: Morris

```bash
cd ~
sbatch $S/sa_run.sbatch design.py morris          # warten, bis der Job fertig ist (squeue)
bash $S/submit.sh results/design_morris.csv       # 2300 Läufe, ≈ 3800 CPU-h (~1.7 h pro Lauf inkl. Auswertung)
```
- `MAXPAR=200 bash $S/submit.sh …` erlaubt mehr gleichzeitige Läufe (Standard 100).
- Mehr Zeit oder Speicher pro Lauf: `SIMTIME=12:00:00 SIMMEM=8G bash $S/submit.sh results/design_morris.csv`.
- Danach laufen automatisch: Sammeln, `evaluate_pooled.py` (gepoolte Ziele pro Design-Punkt), `analyze_morris.py`.

---

## 4. Stufe 2 auswerten und Parameter wählen

Ausgabe `slurmOut/sa_collect-<job>.out`, Dateien in `results/`:

| Datei | Inhalt |
|---|---|
| `morris_indices.csv` | μ*, CI, μ, σ, Rauschschwelle, `in_noise`, `mu_star_prior` pro Zeit × Output × Faktor (`kind`: per run / pooled) |
| `plots/morris_ranking_t*.png`, `morris_musigma_t*.png`, `morris_overview.png` | Rangfolge mit Rauschschwelle, Nichtlinearität, Übersicht |
| `morris_calibration_suggestion.csv` | **Vorschlag** frei / fix pro Parameter (über Rauschen, ≥ 10 % des größten μ*, prior-skaliert) |
| Konsole: Anastomose-Gruppe | μ* von ana, ana_L_scale, R_a, B_a relativ zu distTH |
| `morris_target_correlation.csv`, `plots/target_correlation.png` | Ziele mit |ρ| ≥ 0.9 |
| `morris_state_correlation.csv` | Ziele vs. Wartezeit, Kreuzungen, Zufluss, eingefrorene Spitzen |
| `morris_failure_map.csv`, `plots/morris_failure_map.png` | Ausfälle (excluded / failed) pro Faktor und Stufe |

Kontrollfaktor `R_b_ctrl` muss μ* = 0 haben (sonst stimmt der b_prob-Schalter nicht).

**Entscheiden und eintragen (lokal, dann committen und pullen):**
- In `config.PARAMETERS` bei jedem Parameter `calibrate=True` oder `False` setzen.
- Dominiert distTH die Anastomose-Gruppe: `ana_L_scale`, `R_a`, `B_a` (evtl. `ana`) auf `False`.
- Stark korrelierte Ziele: in `calibrate.TARGETS` eines auskommentieren oder Toleranzen mit `widen` erweitern.
- `R_tropismN` wirkungslos auf `persistence_length`, aber Persistenz weit unter 10 cm: `R_tropismS` als Parameter ergänzen.
- Ausfälle an einem Bereichsende: Bereich begründet enger oder `MAX_WAIT_DAYS` größer.

Prüfen:
```bash
sbatch $S/sa_run.sbatch calibrate.py --check      # "Fixed at the nominal value: ..." zeigt die fixierten Parameter
```

---

## 5. Stufe 3: Kalibrierung

```bash
cd ~ && sbatch $S/calibrate.sbatch --full
```
- Ein Job mit 128 Kernen, 16 Seeds pro Punkt, Zeitlimit 72 h. Pro Iteration 1 Punkt + (Anzahl freier Parameter)
  Differenzpunkte; bei ~1.6 h pro Lauf (120 h nach t*) etwa 5 h pro Iteration, also ~14 Iterationen pro Job.
  `--full` läuft bis zum Minimum statt am ersten Punkt innerhalb der Toleranz zu stoppen.
- Abgebrochen (Zeitlimit, Knotenausfall)? **Einfach dasselbe noch einmal einreichen** – fertige Läufe werden wiederverwendet,
  der Löser spielt deterministisch nach.
- Andere Seeds/Iterationen: `sbatch $S/calibrate.sbatch --full --seeds 24 --max-nfev 40`.

**Ergebnis** in `calibration/`:
| Datei | Inhalt |
|---|---|
| `report.txt` | Status, bester Punkt, Ziele mit F, CHECKS (Ausfälle, Rauschen, Konflikt mit Literatur, ATTENTION Ringe), wodurch jedes Ziel bestimmt ist, Rang der Jacobi-Matrix, **Datenanteil und ±1σ pro Parameter** |
| `evals.csv` | jede Auswertung mit Werten, SE, F |
| `jacobian.csv` | dF/du der letzten Jacobi-Matrix |
| `best_params.xml` | bester Parametersatz |

Lesen: Ziele mit `ok` / `within noise` / `off`; Parameter „<- AT BOUND“ → Bereich zu eng; „determined by the prior“ →
Ziele bestimmen ihn nicht; „Conflict with literature“ → Ziel, Modellstruktur oder Prior prüfen.

---

## 6. Nach der Kalibrierung

- **E7:** `jacobian.csv` mit `morris_indices.csv` vergleichen (lokale vs. globale Rangfolge).
- **Validierung:** Läufe mit `best_params.xml`, Ê_g (`global_efficiency_rel`) und Ĉ (`cost_rel`) bei 20/40/60 h mit
  `config.VALIDATION_OUTPUTS` (0.33/0.44/0.53 und ≈ 0.1) vergleichen. Dafür gibt es noch kein Skript.
- **Sobol** auf den freien Parametern (Ersatzmodell aus Kalibrierungs- und Morris-Läufen) – noch nicht gebaut.

---

## 7. Nützliches

**Nachreichen:** Einfach denselben Befehl noch einmal; fertige Läufe (ok, excluded) werden übersprungen.
```bash
SIMTIME=16:00:00 SIMMEM=8G bash $S/submit.sh results/design_phase0.csv   # fehlende Läufe mit mehr Zeit/Speicher
bash $S/submit_all.sh --no-smoke --no-morris --no-calibration           # ganze Stufe 1 ohne Smoke-Test nachreichen
```

**Nur neu auswerten** (neue Outputs, geänderte Auswertung, ohne Simulation):
```bash
bash $S/submit.sh results/design_morris.csv --eval-only
```

**Einzelne Läufe ansehen:** `results/status_<design>.csv` (ok / excluded / failed / missing),
`runs/<run_id>/meta.json` (Fehlermeldung, Laufzeit, bei excluded `crossed_cm`).

**Status-Begriffe:**
- `ok`: t* erreicht, ausgewertet
- `excluded`: kein t* und ≤ 1 cm Hyphen auf der Pilzseite (Ausschlusskriterium des Experiments), wird nicht wiederholt
- `failed`: kein t* trotz Wachstum (langsam) oder Fehler
- `missing` / `not_evaluated`: Lauf fehlt (Zeitlimit) bzw. nicht ausgewertet

**Jobs hängen mit „DependencyNeverSatisfied“:** Eine Stufe davor ist gescheitert (z. B. Gate) → deren Ausgabe lesen,
wartende Jobs mit `scancel` entfernen, nach der Korrektur neu einreichen.
