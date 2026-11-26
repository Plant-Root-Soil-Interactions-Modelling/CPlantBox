# Sensitivitätsanalyse und Kalibrierung – vollständige Beschreibung der Überarbeitung

#amf #kalibrierung #sensitivitaetsanalyse

Diese Datei beschreibt jede Änderung an SA und Kalibrierung des Hyphenmodells (Petrischale nach Oyarte Galvez et al. 2025,
Nominalmodell `tomatoparameters/TwoHyphaePlusBAS.xml`): das jeweilige Problem, die getroffene Entscheidung, die Umsetzung
im Code und wie man das Ergebnis liest. Am Ende stehen die offenen Schritte und eine Methodenbeschreibung, aus der sich der
Methodenteil schreiben lässt. Die technische Bedienung steht in `README.md` (englisch).

> [!important] Folge für die Läufe
> Das Nominalmodell hat sich geändert (abgeleitete Mindestpfadlänge der Anastomose, `ana` = 0.7, t* über die Runner-Länge,
> Ringe um den Eintrittspunkt bis 32 mm, gemeinsame Parameterliste). Alle bisherigen Läufe sind damit nicht mehr
> vergleichbar: **Phase 0, Langzeittest, Morris, Inokulum-Morris und Kalibrierung müssen neu laufen.** Neue Designs haben
> neue Lauf-IDs (`p0v5_*`, `m4_*`, `inoc5_*`) und die Spalte `factor_set = v3`, damit alte Läufe mit gleicher ID nicht
> stillschweigend wiederverwendet werden. Die erste Morris-Runde bleibt mit ihrer alten Faktorliste (`FACTORS_V1`)
> auswertbar.

> [!important] Simulierte Zeit 120 h nach t*
> Der Langzeittest aus Phase 0 hat gezeigt, dass K1 und K2 erst nach etwa 120 h die Werte des Experiments erreichen
> (nominal K1 763 → 1422 µm mm⁻², K2 0.071 → 0.151 mm⁻² zwischen 60 und 120 h; Ziele 1400 bzw. 0.15). Deshalb gilt
> jetzt `HOURS_HYPHAE = 120`: Die Wellengrößen (alpha, beta, K1, K2, Wellen- und Spitzengeschwindigkeit) werden über
> 0–120 h gefittet und bei 120 h gespeichert, Snapshots bei 10, 20, 30, 40, 60, 90, 120 h. Die Zeitverläufe des Bildgraphen
> bleiben bei 20/40/60 h (Zeitpunkte der Abbildungen). Der Langzeittest läuft jetzt 180 h (`LONG_HOURS`, Bildgraph
> zusätzlich bei 150/180 h) und prüft, ob die Wellengrößen bei 120 h konvergiert sind. Ein Lauf dauert nominal ~1.6 h
> (vorher ~0.4 h): `sim.sbatch` 8 h, Langzeittest 24 h. **Phase 0 muss mit 120 h neu laufen** (Lauf-IDs `p0v5_*`, seit `ana_L_scale` = 3 und Schichtdicke 0.02 cm, siehe A4/A5/A11).

---

## Offene Schritte

### Nach Phase 0
- [ ] `results/phase0_warnings.txt` lesen: Zeilen „ATTENTION …“ zu wenigen Barriere-Kreuzungen bis t* (unsicherer Eintrittspunkt) und zu wachsenden Spitzen jenseits von 32 mm (Ringe zu klein).
- [ ] Konvergenz dx (`phase0_dx.csv`): Weichen die Ziele bei 0.01 cm um mehr als ~2 % von der feinsten Stufe ab? Dann dx verkleinern.
- [ ] dt (`phase0_dt.csv`): dt = 1 h ist festgelegt (A5). Rest-Abhängigkeit am Nominalpunkt mit `ana_L_scale` = 3 dokumentieren (erwartet: Netz bei 0.25 h ~25 % kleiner, Anastomosen pro Spitze gleich).
- [ ] Schichtdicke (`phase0_layer.csv`) und Mindestpfadlänge der Anastomose (`phase0_anaL.csv`): Wie stark reagieren `alpha`, Anastomosen, `encounter_C`? Entscheidung, ob die Annahmen tragbar sind.
- [ ] Common Random Numbers (`phase0_crn.csv`): Verhältnis deutlich unter 1? Sonst brauchen Kalibrierung und Morris mehr Seeds.
- [ ] Streuung zwischen Seeds für die Ziele aus einer Kolonie (Ausgabe „Spread between seeds …“) als `between_cv` bei K1, K2, `front_velocity`, `front_tip_speed`, `alpha`, `beta` in `calibrate.TARGETS` eintragen.
- [ ] Netzgrößen plausibel? Ê_r (`root_efficiency_rel`), Ê_g (`global_efficiency_rel`) zwischen 0 und 1, Ĉ (`cost_rel`) > 0. Sonst ist das Netz auf der Pilzseite nicht zusammenhängend → nur die größte Komponente auswerten.
- [ ] BAS-Anteil: Vergleich Subtyp vs. Bildkriterien (`phase0_bas_share.csv`).
- [ ] Spitzengeschwindigkeit an der Front: verrauscht? Dann `SPEED_FROM` erhöhen.
- [ ] Verzweigungswinkel und Persistenzlänge nominal vs. Experiment (Modell vorab: 52°/57° statt 73°/80°, 4.7 cm statt ≥ 10 cm).
- [ ] `wave_width` stabil definiert? Dann als Ziel aufnehmen (≈ 15 h ± 30 %).

### Nach dem Langzeittest
- [ ] Zufluss neuer Hyphen mit Wurzeln an vs. aus (`long_influx.png`, `long_roots_compare.csv`): Unterscheiden sich die Outputs innerhalb der 120 h? Wenn ja: Wurzeln weiterwachsen lassen.
- [ ] Wellengrößen bei 120 h vs. 150/180 h (`long_summary.csv`): konvergiert? Sonst `HOURS_HYPHAE` weiter erhöhen.

### Nach der Morris-Runde
- [ ] Vorschlagstabelle (`morris_calibration_suggestion.csv`) lesen und in `config.PARAMETERS` von Hand `calibrate=True/False` setzen.
- [ ] Anastomose-Gruppe (distTH, ana, ana_L_scale, R_a, B_a): Dominiert distTH, die anderen fixieren.
- [ ] Zielkorrelation (`morris_target_correlation.csv`): Gruppen mit |ρ| ≥ 0.9 → ein Ziel genügt oder Toleranzen um √(Anzahl) erweitern.
- [ ] Ausfallkarte (`morris_failure_map.png`): Häufen sich Ausfälle an einem Ende eines Bereichs, Bereich begründet enger oder `MAX_WAIT_DAYS` größer.
- [ ] Zustand bei t* (`morris_state_correlation.csv`): Welche Ziele hängen an der Wartezeit? Welche Parameter wirken über die Wartezeit?
- [ ] Inokulum-Morris (`inoculum_vs_hyphae.csv`): Verhältnis ≪ 1 → Normierung auf t* nimmt die Wurzelphase heraus.
- [ ] Bei `R_tropismN` an der Obergrenze: `R_tropismS` als Parameter aufnehmen.

### Nach der Kalibrierung
- [ ] E7: Jacobi-Matrix (`jacobian.csv`) mit Morris vergleichen (lokale vs. globale Rangfolge, Nichtlinearität).
- [ ] Validierung: Läufe mit `best_params.xml`, Ê_g und Ĉ mit `config.VALIDATION_OUTPUTS` vergleichen (Skript noch zu bauen).
- [ ] Varianzbasierte SA (Sobol über Ersatzmodell) auf den Parametern, die Morris als wirksam zeigt.

### Unabhängig davon
- [x] Definition der Anastomosendichte ρ_A im Paper: **keine Angabe** (Einheit, Normierung) → kein Ziel; ohnehin redundant zu Anastomosen und L_RH.
- [ ] Absterben: Kinder toter Hyphen werden im C++ nicht weiter simuliert (A6) – Korrektur bewusst zurückgestellt.

---

## A – Annahmen im Modell-Setup

### A1 Wurzelwachstum wird bei t* abgeschaltet
**Problem:** Ab t* wachsen die Wurzeln nicht mehr, Kolonisierung und Hyphenentstehung laufen weiter. Neue Hyphen kommen nur
noch von den vorhandenen Wurzeln; im Experiment wachsen die Wurzeln weiter. Die Vermutung: irgendwann kommen im Modell keine
neuen Hyphen mehr über die Barriere. Ob das innerhalb der ausgewerteten Zeit eine Rolle spielt, war unbekannt.

**Lösung:**
- Diagnose des Zuflusses (`model_hook.influx_diagnostics`), pro Stunde über die ganze Simulation aus der fertigen Pflanze:
  neue Hyphen aus Wurzeln, Segmente, die die Barriere nach innen bzw. außen überqueren (gesamt, Runner, BAS), neue Wurzellänge.
- Outputs `influx_tstar` (Kreuzungen pro Stunde in den 6 h vor t*), `influx_window` (Mittel 0–120 h), `influx_decay`
  (letzte 10 h / bei t*: 1 = gleichbleibender Zufluss, 0 = versiegt).
- Langzeittest (`design.py long`): 180 h nach t* (`LONG_HOURS`), 10 Seeds mit Wurzeln aus und dieselben Seeds mit weiterwachsenden Wurzeln
  (`roots_on`). `analyze_long.py` vergleicht die Outputs gepaart über die Seeds bei 60/90/120 h und plottet den Zufluss
  (`plots/long_influx.png`); es gibt aus, ab wann der Zufluss unter 10 % seines Maximums fällt.

**Lesen:** Unterscheiden sich die Outputs innerhalb des ausgewerteten Fensters nicht, ist die Annahme für SA und
Kalibrierung unkritisch, unabhängig davon, was danach passiert.

**Ergebnis (Phase 0, Lauf über 120 h):** Wurzeln an vs. aus ohne signifikanten Unterschied – Annahme bleibt. Dabei hat
sich gezeigt, dass K1 und K2 erst nach ~120 h die Zielwerte erreichen → `HOURS_HYPHAE = 120`, Langzeittest 180 h.

### A2 Wurzelphase, Inokulum, Wurzelarchitektur
**Problem:** Dauer der Wurzelphase, Lage des Samens und die Wurzelarchitektur sind nicht belegt. Die Wurzelarchitektur ist
eine modifizierte Tomate als Ersatz für die stark verzweigte Wurzelorgankultur (ROC) von *D. carota* des Experiments, ohne
Daten zur Geometrie. Diese Größen bestimmen, wann und wo Hyphen durch die Öffnung kommen – die Randbedingung des ganzen
Pilzkompartiments.

**Entscheidung:** Entstehungsdichte (4/cm), `lmbd` (0.15), `vi`, `hyphalDelay`, `maxColonization`, `minAge`, `maxAge` stammen aus
Schnepf et al. und bleiben fest (Quellen jetzt im Code und im XML). Unbelegt sind die Wurzelphase (`ROOT_DAYS` = 5 d: „die
Wurzeln sollen erst etwas wachsen“), die Lage des Wurzelstücks und die Architektur.

**Lösung:** Eigene Inokulum-Morris (`design.py inoculum`, `config.INOCULUM_FACTORS`), Hyphenparameter nominal:
| Faktor | Bereich | Wirkt auf |
|---|---|---|
| `root_days` | 3–10 d | Dauer der Wurzelphase |
| `root_r_scale` | ×/÷ 2 (log) | Wachstumsrate `r` aller Wurzeln |
| `root_ln_scale` | ×/÷ 2 (log) | Verzweigungsabstand `ln` aller Wurzeln |
| `seed_pos_x` | 5 mm vom Rand bis ~1 cm vor der Barriere | Lage des Wurzelstücks |
| `emergence_density` | 2–6 /cm (Schnepf ±50 %) | Vergleich mit einer belegten Größe |

10 Trajektorien × 6 Punkte × 3 Replikate = 180 Läufe. Outputs zusätzlich: `wait_days`, Zufluss, `frozen_share`.
`analyze_morris.py inoculum` schreibt `inoculum_vs_hyphae.csv`: größtes μ* der Inokulum-Faktoren / größtes μ* der
Hyphenfaktoren pro Output.

**Lesen:** Verhältnis ≪ 1 → die Normierung auf t* nimmt den Einfluss der Wurzelphase heraus, die Tomate als Ersatz ist für
die Ergebnisse unkritisch. Bekannte Einschränkung: Eine ROC wird als Wurzelstück angesetzt, die Hauptwurzel der Tomate ist ein
Artefakt.

### A3 Definition von t*
**Entscheidung:** 10 cm Hyphenlänge hinter der Barriere wie Oyarte Galvez et al. 2025 (gewählt, um Ungenauigkeiten bei der
visuellen Bestimmung des Übertritts zu vermeiden). Im Paper ist t* erreicht, wenn die **Runner-Länge** L_RH = 10² mm ist.

**Lösung:** `makeSimulation(..., crossing_subtype = 1)` zählt nur die Runner-Länge (vorher alle Hyphen inkl. BAS, wodurch t*
zu früh kam). Bleibender Unterschied: Modell 3D-Länge, Experiment projizierte Länge im Bild.

### A4 Dicke der Hyphenschicht
**Problem:** Die Hyphen wuchsen im Modell in einer 1 mm dicken Schicht (3D), auf dem Cellophan quasi in 2D. Bei
`distTH` ~100 µm bestimmt die Schichtdicke die Begegnungswahrscheinlichkeit und damit Anastomosen, `alpha`, `encounter_C`.

**Befunde:**
- Phase 0 (p0v2, L×1, 60 h): 0.1 → 0.02 cm gibt `alpha` −65 %, `encounter_C` −63 %, Gesamtlänge +64 % (dickere Schicht:
  weniger Begegnungen, größeres Netz).
- Lokaler Schichttest (2026-11, Nominalmodell mit `ana_L_scale` = 3, bis t* + 48 h, 3–7 Seeds je Schicht, gepaart; nach
  23 von 32 Läufen abgebrochen). Verhältnis zu 0.005 cm:
  - 0.02 cm: Anastomosen pro Spitze 1.00 [0.83, 1.20], Verzweigung 0.99, Netzlänge 0.94
  - 0.01 cm: Anastomosen pro Spitze 0.81, Netzlänge 0.80
  - 0.1 cm: Anastomosen pro Spitze 0.61 [0.53, 0.70], Netzlänge 1.51 [1.27, 1.80]

  Ab 0.02 cm abwärts ändert sich also nichts Systematisches mehr (quasi-2D-Grenze), 0.1 cm liegt deutlich daneben. Die
  Laufzeit hängt nicht an der Schichtdicke.
- Die Schicht hält: Auf der Pilzseite (x > b/2) liegen keine Spitzen oder Knoten unterhalb der Schicht oder außerhalb der
  Geometrie.
- Die Öffnung in der Barriere ist im Modell 2 mm tief, die Schicht dahinter nur h. Spitzen, die tief in die Öffnung
  wachsen, laufen bei x = b/2 gegen die Stufe und bleiben dort (bei 0.02 cm etwa ein Drittel der aktiven Spitzen nahe
  der Barriere, davon 71 % unterhalb der Schicht). Sie liegen außerhalb der Ringe und zählen damit nicht in den
  Ring-Outputs.

**Entscheidung:**
- Nominal `HYPHAE_LAYER` = 0.02 cm.
- Die Schichtdicke ist ein Faktor der SA: `hyphae_layer` in `config.SA_GEOMETRY`, 0.005–0.1 cm, log, nur in der SA,
  nicht kalibriert. Er geht als Argument an die Simulation (`target="sim"`).
- Die Öffnung bleibt, wie sie ist. Eine Rampe oder eine flachere Öffnung entspräche nicht dem Versuchsaufbau.
- Die steckengebliebenen Spitzen werden nur beobachtet: `raw.npz` enthält `diag_opening_tips`,
  `diag_opening_tips_below_layer` und `diag_fungal_tips_end`. Ein Effekt der Schichtdicke über die Wartezeit bis t* und
  den Zustrom bleibt Teil ihres Morris-Effekts.
- Phase-0-Gruppe: 0.02 (nominal) / 0.005 / 0.1 cm, also Nominalwert und die Enden des Morris-Bereichs
  (`phase0_layer.csv`).

### A5 Zeitschritt dt = 1 h
**Problem:** 1 h ist das Bildintervall des Experiments, also das Auswerteintervall, nicht zwingend der numerische Schritt.
Die Anastomose wird einmal pro Schritt am Spitzenknoten geprüft (`MycorrhizalPlant.cpp:304`); bei 240 µm/h springt die
Spitze pro Schritt 240 µm > `distTH`. Auch `ana` (pro Schritt) und die BAS-Teilung (b·dt pro Schritt) hängen an dt.
Phase 0 (v2) zeigte: Bei 0.25 h ist das Netz 6–7× kleiner.

**Lösung:**
- Konvergenztest in Phase 0: 24 / 48 / 96 Schritte pro Tag (1 / 0.5 / 0.25 h), Auswertung weiter stündlich (`phase0_dt.csv`).
- Fehler behoben: Die Hyphenphase lief über `hours_hyphae` *Schritte*; jetzt über Stunden (`round(hours / (dt · 24))`).

**Isolierte Tests (November 2026, 8 Seeds):** Gleicher Zustand bei dt = 1 h (Wurzelphase + 2 d Hyphen), dann Wurzeln und
Kolonisierung aus und 2 d Hyphen bei 1 / 0.5 / 0.25 h; gemessen pro aktiver Spitze und Tag.
- Verzweigung (Runner-Seitenäste im Abstand ln, BAS-Spitzenteilung) und Wurzelphase hängen nicht von dt ab.
- Anastomose: Mit der abgeleiteten Mindestpfadlänge (×1) verschmelzen junge Seitenäste (~1 h alt, ~190 µm lang) mit der
  eigenen Eltern- bzw. Kindhyphe nahe dem Verzweigungspunkt: Sie biegen durch den Tropismus parallel zur Elternhyphe
  zurück (wirksamer Winkel ~6° statt θ) und bleiben innerhalb von distTH, sobald der Pfad länger als L ist. Anteil an
  allen Anastomosen 38 % bei 1 h, 91 % bei 0.25 h; Rate pro Spitze 0.32 → 2.67 /d. Bei 1 h überspringt die Spitze den
  Bereich meist (240 µm pro Schritt, SDF einen Schritt alt), daher wirkte dt = 1 h wie ein längeres L.
- Ab 3 × L ist die Anastomose-Rate pro Spitze unabhängig von dt (0.25 / 0.23 / 0.22 /d), Verwandte nur noch ~0.05 /d
  (echte Schleifen im eigenen Baum).
- Rest: Ohne Anastomose ist das Netz nach 2 d bei 0.25 h ~15 % kleiner (Wachstum pro Spitze −6 %, signifikant; die
  Seitenäste folgen dem Wachstum), mit Anastomose ~25 % (vermutlich Weiterwachsen nach dem Kontakt bis Schrittende und
  ~15 % übersprungene Kontakte bei 1 h).

**Entscheidung:** dt = 1 h gehört zur Modelldefinition (Bildintervall des Experiments); die kalibrierten Parameter gelten
für dt = 1 h. Das Artefakt der Selbst-Anastomose ist über `ana_L_scale` = 3 behoben (A11). Die Rest-Abhängigkeit des
Wachstums von dt ist dokumentiert, nicht korrigiert (Ursache vermutlich Zuwachs pro Schritt < dx oder Geometrie der
dünnen Schicht). Eine C++-Korrektur (Anastomose einmal pro Begegnung, ganzer Weg geprüft) bleibt möglich.

### A6 Artefakt beim Absterben
**Problem:** `Hyphae.cpp:108` simuliert die Kinder einer Hyphe nur, solange sie lebt. Stirbt sie nach `hlt`, friert der
ganze Teilbaum ein, auch jenseits der Barriere. Die Sensitivität auf `hlt` mischt damit Biologie und Artefakt; bei langer
Wartezeit vor t* betrifft das schon die Kolonie bei t*.

**Lösung (Diagnose, keine Code-Korrektur):** `model_hook.frozen_diagnostics` bestimmt pro Stunde abgestorbene Hyphen und
eingefrorene Spitzen (gesamt und auf der Pilzseite); Output `frozen_share` (Anteil eingefrorener an allen Spitzen auf der
Pilzseite am Ende). Am echten Modell geprüft: Die Diagnose findet genau die Spitzen, die eine unabhängige Zählung findet.

**Entscheidung (November 2026, `frozen_share` nominal 0.41):** bleibt so. Biologisch wachsen Kinder einer toten Hyphe nicht
weiter, weil sie nicht mehr versorgt werden. Wie viel davon betroffen ist, ist eine Frage der Definition von `hlt` und
`lmax` (Lebensdauer und Endlänge), nicht der Numerik. Root macht es genauso (`Root.cpp:105`, Kinder nur innerhalb von
`if (alive)`), dort wirkt es nur nicht, weil `rlt` ≈ 1e11 d ist.

**Mögliche Erweiterung (nicht umgesetzt):** ein echtes Absterben durch die Elternhyphe statt Einfrieren. In
`Hyphae::simulate` bei `alive = false` (Lebensdauer erreicht) rekursiv alle Nachfahren auf `alive = false` und
`active = false` setzen, eventuell mit einer Verzögerung (Restversorgung) als Parameter. Dazu beim eigenen Tod `active =
false` setzen: Heute behält eine tote oder eingefrorene Spitze `active = true` und wird in
`MycorrhizalPlant::simulateAnastomosis` (`h->isActive()`) weiter auf Anastomose geprüft, obwohl sie nicht mehr wächst.
Die Auswertung zählt wachsende Spitzen schon über die Knotenzeiten (`AMFDensities.getActiveTipsPerRingPerSubtype`), nicht
über `isActive()`.

### A7 Radius der Runner
**Entscheidung:** 50 µm **Radius** (Schnepf et al. 2017, dort nur Runner). Literatur-Prior 25–100 µm (Faktor 2, ±2σ);
Untergrenze des Suchbereichs 10 µm (Bago et al. 1998: Durchmesser > 20 µm). Hinweis: Im Modell wirkt `a` über
„Abstand minus Radius“ direkt auf die Anastomose, anders als bei Schnepf.

### A8 Struktur der BAS
**Problem:** Ein BAS-Segment teilt sich nur, solange es wächst (bis `lmax`); nominal ~26 % Teilungswahrscheinlichkeit pro
Segment, also < 1 Kind (gewollt: BAS leben nicht lang). Die Größe der BAS hängt damit an b · lmax / v; `lmax` war fest
(„≈ 3·dx“, rein numerisch begründet).

**Lösung:** `B_lmax` (0.01–0.1 cm, log, schwacher Prior) wird kalibriert und in der SA variiert. `la`, `lb`, Tropismus der
BAS bleiben fest.

### A9 Streuungen (dev)
**Problem:** Die Streuungen der Parameter sind Schätzwerte ohne Bereiche.

**Lösung:** SA-Faktor `cv_scale` (×/÷ 2, log) skaliert alle `dev` außer dem von `v` gemeinsam. Ohne Wirkung → die
Schätzungen sind unkritisch.

### A10 Kleinere Punkte
- **Ringauflösung:** `analyze_rings.py` fasst je 3 Ringe zusammen (15 → 5, wieder gleiche Fläche) und vergleicht alle
  ringabhängigen Outputs pro Lauf und gepoolt, ohne Neusimulation.
- **root_dist** = 5 mm wie Oyarte Galvez (Kommentar ergänzt).
- **Anastomose im eigenen Baum:** Das Modell erlaubt Fusion im eigenen Baum ab einer Mindestpfadlänge entlang des Myzels.
  `AMFDensities.getTipDistances` schloss den ganzen eigenen Baum aus; jetzt wie im Modell über die Pfadlänge
  (`exclude="path"`, alte Logik mit `exclude="tree"`). Am echten Modell: 786 statt 40 von 855 Spitzen mit Kandidat
  innerhalb 1 mm – die alte Version hat Begegnungen stark unterschätzt.

### A11 Mindestpfadlänge der Anastomose (`anaMinDist`)
**Problem:** Nicht im XML gesetzt, also galt der Default 0.1 cm aus `hyphaeparameter.h`; gedacht war die abgeleitete Formel.
Die abgeleitete Formel distTH / sin θ + 2·dx (≈ 0.032 cm) nimmt an, dass ein Seitenast gerade unter θ wegwächst. Das tut
er nicht (A5): Er biegt parallel zur Elternhyphe zurück und verschmilzt mit ihr.

**Lösung:** Faktor `ana_L_scale` auf den abgeleiteten Wert in SA und Kalibrierung, **nominal 3** (≈ 0.095 cm, fast der alte
Default), Bereich und Prior 1.5–6 (log, ×/÷ 2). Im Template-XML steht der Wert explizit (0.09464 cm), `write_xml` setzt
ihn aus dem Faktor neu (skaliert mit distTH, θ, dx). Phase-0-Gruppe ×3 (Referenz) / 1 / 1.5 / 6, `phase0_anaL.csv`.

### Ringe (aus F6)
**Problem:** Oyarte Galvez berechnen die radialen Dichten über den Abstand r vom **Barriere-Eintrittspunkt**, pro Replikat.
Das Modell hatte 15 Ringe gleicher Fläche über die ganze Schale (47 mm) um die Schalenmitte. Die Punkte in Fig. 2c
(r ≈ 11–31 mm) passen zu 15 Ringen gleicher Fläche bis ≈ 32 mm; die Modellringe waren ~1.5× breiter. Ringbreite verlängert
die Wellendauer und verflacht K2.

**Lösung:**
- `AMFSimulationSetup.crossingPoint`: Median der Schnittpunkte aller nach innen kreuzenden Segmente bis t* = Eintrittspunkt.
- `AMFSimulationSetup.ringsAround`: 15 Ringe gleicher Fläche bis `RING_RADIUS` = 3.2 cm um diesen Punkt, Flächen numerisch
  (bei verschobenem Zentrum schneidet der Schalenrand die äußeren Ringe).
- `model_hook` baut die Ringe pro Lauf nach der Simulation; alle Ringgrößen (K1, K2, Wellenfits, alpha/beta, Frontposition,
  „Front“ der Spitzengeschwindigkeit) nutzen sie. `RING_CENTER = "dish"` stellt die alten Ringe wieder her.
- Warnungen: weniger als 5 Kreuzungen bis t* (Median) oder > 10 % der wachsenden Spitzen jenseits 32 mm →
  „ATTENTION …“ in `results/phase0_warnings.txt` und `calibration/report.txt`.

---

## B – Verhältnis von SA und Kalibrierung

### B11 Kalibrierziele in der Morris-Analyse
**Problem:** Die Morris-Outputs enthielten fast keines der Kalibrierziele; pro Lauf scheitern die Wellenfits zu oft.

**Lösung:** `evaluate_pooled.py` poolt pro Design-Punkt über die Replikate wie die Kalibrierung über ihre Seeds (`alpha`,
`beta`, K1, K2, `front_velocity`, `wave_width`), Standardfehler per Bootstrap (100). Läuft automatisch in `collect.sbatch`.
Im Dummy-Test an 265 von 280 Punkten definiert.

### B12 Einheitliche Definitionen
**Lösung:** Eine Definition für SA, Kalibrierung und alle Auswertungen (`evaluate.prepare`, `evaluate.wave_quantities`):
- Spitzen = wachsende **Runner**-Spitzen (wie im Experiment).
- `front_velocity` aus der Ankunftszeit t_n des Sigmoid-Fits der Fadendichte (Theil-Sen); Variante aus dem Spitzenpeak heißt `v_tips`.
- Nomenklatur: **alpha = Anastomose** (a = α ρ n, ~23 µm/h), **beta = Verzweigung** (b = β n, ~4 %/h); im Paper umgekehrt.
  `AMFDensities.anastomosisAlpha` (alter Name `anastomosisBeta` bleibt), Plot-Beschriftungen angepasst.

### B13 Ein gemeinsamer Parameterraum
**Lösung:** `config.PARAMETERS` definiert jeden Parameter einmal: Bereich (= Morris-Bereich = Suchbereich), Prior, Quelle,
`calibrate` (frei oder fest). Die SA variiert zusätzlich `config.SA_ONLY` (`cv_scale`, Kontrollfaktor `R_b_ctrl`).
`xml_tools.write_xml` kann `successor`, die Kopplung lb = ln, `v_cv`, `ana_L_scale` und `cv_scale` schreiben.

| Parameter | Bereich | Prior / Quelle |
|---|---|---|
| R_v, B_v | 0.432–0.72 cm/d | 180–300 µm/h (Oyarte Galvez) |
| R_hlt | 4–12.5 d | 4–12.5 d (Schnepf 2016/2008) |
| R_lmax | 1.5–5 cm | schwach |
| R_ln (= lb) | 0.05–0.6 cm, log | 0.1–0.6 cm |
| R_theta, B_theta | 40–100° | schwach (Winkel sind Ziele) |
| R_tropismN | 1–8 | schwach (Persistenzlänge ist Ziel) |
| R_succ_BAS | 0.5–0.95 | schwach |
| R_a | 10–100 µm Radius, log | 25–100 µm (Schnepf 2017) |
| B_b | 0.3–2 /d, log | schwach (Messung ist Ziel beta) |
| B_hlt | 3–12 d | 5–9 d (Bago 1998, Kokkoris 2026) |
| B_tropismN | 0.5–4 | schwach |
| B_a | 0.5–50 µm Radius, log | 1–5 µm (Bago 1998) |
| B_lmax | 0.01–0.1 cm, log | schwach |
| ana | 0.3–1 (nominal 0.7) | schwach |
| distTH | 30–500 µm, log | 50–300 µm (Richter et al. 2024) |
| ana_L_scale | 1.5–6, log (nominal 3) | 1.5–6 |
| v_cv | 0.15–0.40 | 0.2–0.4 (gemessener CV ≈ 0.29) |
| cv_scale (nur SA) | 0.5–2, log | – |
| R_b_ctrl (nur SA, Kontrolle) | – | soll μ* = 0 geben |

### B14 Reihenfolge
1. Phase 0 und Langzeittest → Numerik und Strukturannahmen festlegen.
2. Morris-Screening (Hyphen) und Inokulum-Morris.
3. Parameter **von Hand** auswählen: `analyze_morris.py` schreibt nur einen Vorschlag.
4. Kalibrierung der freien Parameter, übrige auf Nominalwert (stehen in `report.txt`).
5. Unsicherheit und varianzbasierte SA um den kalibrierten Punkt.

---

## C – Kalibrierung

### C15 Herkunft und Definition der Ziele
- **BAS-Anteil:** Das Experiment klassifiziert Kanten des Bildgraphen als BAS (Kante < 400 µm, oder < 1000 µm mit einer
  Spitze als Endpunkt, oder Breite < 7 µm und Breite × Länge < 9000 µm²). `AMFNetwork.basShareImage` wendet dieselben
  Kriterien auf den ebenen Modellgraphen an (Breite = 2·Radius des Organs, längengewichtet). Am echten Modell 39 % nach
  Bildkriterien gegenüber 16 % nach Subtyp – Kreuzungen und Seitenäste zerstückeln die Runner in kurze Kanten.
- **beta:** zählt neue Runner-Spitzen (aus Runner-Seitenästen); bestimmt R_ln, R_v, R_succ_BAS, nicht B_b. Nominal ≈ 6 %/h.
- **Welle 280 > Spitzen 240 µm/h:** Die Richtung kann die Front nur verlangsamen; schneller wird sie durch die schnellen
  Spitzen (Streuung) → v_cv wird kalibriert.

### C16 Ziele und Toleranzen
Alle Werte aus den Abbildungen von Oyarte Galvez et al. 2025 abgelesen (vergrößert ausgeschnitten), ungefähr; Quelle steht
pro Ziel im Feld `source`.

| Ziel | Wert | Toleranz | Quelle / Bemerkung |
|---|---|---|---|
| alpha | 23 µm/h | ±30 % ⊕ between_cv | Fig. 3b, ein Replikat; Band von a(t) ±40–50 % pro Zeitpunkt |
| beta | 0.04 /h | ±30 % ⊕ between_cv | Fig. 3b, ein Replikat; Band von b(t) ±30–40 % |
| K1 | 1400 µm mm⁻² | ±12 % ⊕ between_cv | Fig. 2c, eine Kolonie |
| K2 | 0.15 mm⁻² | ±33 % ⊕ between_cv | Fig. 2b, eine Kolonie |
| front_velocity | 280 µm/h | ±10 % ⊕ between_cv | Fig. 2d Inset, ein Fit |
| front_tip_speed | 300 µm/h | ±10 % ⊕ between_cv | Fig. 2d, Spitzen an der Front |
| branch_angle_rh | 73° | ±10 % | Ext. Data Fig. 11b |
| branch_angle_bas | 80° | ±10 % | Ext. Data Fig. 11b |
| persistence_length | ≥ 10 cm | Untergrenze (50 %) | Ext. Data Fig. 11a |
| L_RH_image 20/40/60 h | 280 / 500 / 850 mm | 30/30/40 % ×√3 | 12 Replikate |
| bas_share_image 20/40/60 h | 0.37 / 0.385 / 0.40 | 14/13/13 % ×√3 | 12 Replikate |
| anastomoses 20/40/60 h | 75 / 240 / 445 (seit t*) | 50/60/60 % ×√3 | 12 Replikate, 21 bei t* abgezogen |
| root_efficiency_rel 20/40/60 h | 0.69 / 0.665 / 0.63 | ±0.1 / 0.08 / 0.08 ×√3 | 12 Replikate |

- **Zeitreihen ×√3** (`widen`): Die drei Zeitpunkte eines Replikats sind korreliert; zusammen sollen sie wie ein Ziel zählen.
  Umgesetzt als F / √3 (gleichwertig für logarithmische und lineare Residuen).
- **Eine Kolonie** (`replicates = 1`): Fig. 2 und Fig. 3 zeigen je ein einziges Replikat (betrifft K1, K2,
  front_velocity, front_tip_speed, alpha, beta); die Bänder in Fig. 2 sind räumliche Bootstraps innerhalb
  eines Rings (10 000 Rechtecke), nicht die Streuung zwischen Platten. Wirksame Toleranz √(tol² + between_cv²), between_cv aus
  Phase 0 (CV zwischen Seeds am Nominalpunkt). Bis dahin 0, diese Ziele sind also noch zu streng.
- **Untergrenze** (`bound = "lower"`): Oberhalb von 10 cm ist F = 0.
- **Nicht als Ziel:** L_BAS (folgt aus L_RH und ρ_BAS), ρ_A (im Paper nicht definiert; beide plausiblen Varianten nur als
  Output), Ê_g und Ĉ (zur Validierung zurückgehalten).

### C17 Identifizierbarkeit
`report.txt` zeigt pro Parameter den **Datenanteil** = (JᵀJ)_jj / ((JᵀJ)_jj + (PᵀP)_jj): wie viel Information von den Zielen
statt vom Prior kommt; < 0.5 wird als „durch den Prior bestimmt“ markiert.

### C18 Ergebnis und Unsicherheit
- ±1σ pro Parameter am besten Punkt aus (JᵀJ + PᵀP)⁻¹, linearisiert, in physikalischen Einheiten; σ = Änderung, die die Ziele
  um etwa eine Toleranz verschiebt.
- `--full`: nicht am ersten Punkt innerhalb der Toleranz stoppen, sondern bis zum Minimum laufen (für den finalen Lauf).
- Bekannte Einschränkung: `NAN_PENALTY` und `MIN_OK` machen F unstetig.

### C19 Common Random Numbers
Phase-0-Gruppe: `R_v` um den Kalibrierschritt (0.15 normiert) verschoben, einmal mit den Seeds der Rausch-Läufe, einmal mit
unabhängigen Seeds; `phase0_crn.csv` gibt sd(Differenz, gleiche Seeds) / sd(Differenz, unabhängige Seeds) aus.

---

## D – Methodik der SA

### D20 distTH und ana
- **distTH** ist biologisch begründet größer als der Hyphendurchmesser: Spitzen scheinen andere Hyphen zu „sehen“; ~100 µm
  plausibel, möglicherweise mehr (Richter et al. 2024). Prior 50–300 µm, Bereich 30–500 µm.
- **ana:** unbekannt, Fusion passiert nicht immer → nominal 0.7 statt 1 (vorher am Rand), Bereich 0.3–1.
- **Austauschbarkeit:** distTH, ana, ana_L_scale, R_a, B_a wirken fast nur als Produkt auf die Anastomoserate.
  `analyze_morris.py` gibt ihre μ* für die Anastomose-Ziele relativ zu distTH aus; dominiert distTH, die anderen fixieren.
- **R_hlt:** Suchbereich auf den Prior erweitert (4–12.5 statt 4–10 d).

### D21 Auswertung der Morris-Analyse
- **Schwelle:** statt „2 × Rauschniveau“ das 95 %-Quantil der μ*-Verteilung bei reinem Rauschen, per Monte Carlo aus den
  Standardfehlern der Elementareffekte jeder Trajektorie (Replikate bzw. Bootstrap). Ungleich verteiltes Rauschen geht ein.
- **Bereichsabhängigkeit:** zusätzlich `mu_star_prior` = μ* × (Prior-Breite / Bereichsbreite): Wirkung der tatsächlichen
  Unsicherheit. Vorschlag „fix (prior narrow …)“, wenn ein Parameter über den Bereich wirkt, seine echte Unsicherheit aber
  kaum zählt.
- **Screening:** Morris bleibt Screening (SALib). Wie viele Parameter in den Sobol-Schritt gehen, entscheidet die Morris-Runde.

### D22 Ausfälle
- SALib bleibt (vollständige Trajektorien). Zusätzlich die **Ausfallkarte** (`morris_failure_map.png`,
  `morris_failure_map.csv`): pro Faktor und Morris-Stufe der Anteil der Läufe ohne Ergebnis und der Design-Punkte mit
  undefinierten Zielen; Faktoren mit Häufung an einem Ende werden aufgelistet.
- **Ausschluss wie im Experiment:** Platten mit nicht mehr als 1 cm Hyphenlänge im Hyphenkompartiment wurden als ungesund
  ausgeschlossen. Läufe ohne t* mit ≤ `EXCLUDE_CROSSED_CM` = 1 cm auf der Pilzseite erhalten den Status „excluded“,
  die übrigen „failed“ (langsam). Ausgeschlossene Läufe werden beim erneuten Einreichen nicht wiederholt.

### D23 Zustand der Kolonie bei t*
- `wait_days`, `n_crossings_tstar`, `influx_tstar`, `frozen_share` sind SA-Outputs jeder Morris-Analyse
  (`config.STATE_OUTPUTS`): Welche Hyphenparameter wirken schon vor t*?
- `morris_state_correlation.csv`: Spearman-Korrelation jedes Ziels mit diesen Größen über alle Läufe. Hoch → der Zustand bei
  t* prägt das Ziel, ein Parameter mit großem μ* auf die Wartezeit wirkt dann wahrscheinlich indirekt.
- Keine Daten zur Wartezeit im Experiment → kein Ziel.

---

## E – Auswertungen ohne Neusimulation
- **E1** gepoolte Ziele pro Morris-Punkt ✓ (B11) · **E2** wait_days ✓ (D23) · **E4** Snapshot-Zeiten 10, 20, 30, 40, 60, 90, 120 h ✓ ·
  **E8** Langzeittest ✓ (A1)
- **E6** `morris_target_correlation.csv` und `plots/target_correlation.png`: Spearman-Korrelation der Kalibrierziele über die
  Design-Punkte (Mittel der Replikate bzw. gepoolter Wert), Paare mit |ρ| ≥ 0.9 werden aufgelistet.
- **E3** Diagnosen pro Lauf (keine Ziele): `ana_share_bas` (Anteil der Anastomosen durch BAS-Spitzen), `tip_share_bas`
  (Anteil der BAS an den wachsenden Spitzen).
- **E7** nach der Kalibrierung; **E5** (räumliche Frontfits) weggelassen, redundant.

## F – Zusätzliche Zielgrößen
- **F1/F2 Zeitverläufe und Anastomosen:** Bildgraph zu früheren Zeiten aus der fertigen Pflanze (`t_max` in
  `getPlanarNetwork`, `basShareImage`, `getTerminalRootNodes`, `relativeRootEfficiency`), berechnet bei 20, 40, 60 h
  (`IMAGE_HOURS`) → Ziele in C16.
- **F3 Spitzengeschwindigkeit an der Front:** `model_hook.tip_speed_diagnostics` – pro Stunde und wachsender Runner-Spitze
  gewachsene Bogenlänge / Zeit (exakte Knotenzeiten); Front = äußerste 10 % (Abstand vom Eintrittspunkt); Mittel ab 10 h.
  Am echten Modell Mittel aller Spitzen 245 µm/h (passt zu v = 240).
- **F4 Winkel und Persistenzlänge:** `branching_angle_diagnostics` (Winkel zwischen Anfangsrichtung des Kindes über 200 µm und
  Tangente der Elternhyphe über ±100 µm, projiziert, getrennt RH/BAS) und `persistence_diagnostics` (⟨cos Δθ(ΔL)⟩ bis 1 cm,
  Fit exp(−ΔL/ℓ_p), Median über die Runner).
- **F5 Ê_g und Ĉ:** `relativeRootEfficiency(global_eff=True)`, Normierung (X − X_MST)/(X_DT − X_MST), Ê_g über höchstens 300
  zufällige Knoten aus Spitzen, BAS-Ansätzen und Wurzelknoten; **zur Validierung zurückgehalten** (`config.VALIDATION_OUTPUTS`:
  Ê_g 0.33/0.44/0.53, Ĉ ≈ 0.1).
- **F6 Wellendauer:** Ringgeometrie angepasst (siehe Ringe); als Ziel erst nach Phase 0.

---

## Methode (für den Methodenteil)

### Modell und Versuchsaufbau
- CPlantBox mit Hyphenerweiterung: Runner mit Seitenästen, BAS mit Tip-Splitting; Anastomose, wenn eine Spitze näher als
  distTH an einer Hyphe ist (mit Wahrscheinlichkeit ana), nur ab einer Mindestpfadlänge entlang des Myzels.
- Geteilte Petrischale wie Oyarte Galvez et al. 2025 (94 mm, Barriere 1.6 mm, Öffnung 50 × 2 mm); Hyphen auf der Pilzseite in
  einer 0.2 mm dicken Schicht an der Oberfläche (quasi 2D, Cellophan; in der SA 0.05–1 mm).
- Wurzelorgankultur von *D. carota* durch eine modifizierte Tomatenarchitektur ersetzt; Wurzelphase 5 d, danach Kolonisierung
  und Hyphenentstehung nach Schnepf et al.; Wurzelwachstum ab t* aus.
- t*: 10 cm Runner-Länge hinter der Barriere (wie L_RH = 10² mm); danach 120 h, dt = 1 h (Wellenfits über 0–120 h, Zeitverläufe bei 20/40/60 h).
- Auswertung wie im Experiment:
  - radiale Dichten in 15 Ringen gleicher Fläche bis 32 mm um den Barriere-Eintrittspunkt des jeweiligen Laufs
  - Spitzendichte aus wachsenden Runner-Spitzen
  - Wellenfits: Sigmoid für die Fadendichte, Glocke für die Spitzendichte, Qualitätskriterien wie im Paper
  - ebener Bildgraph (Projektion, Kreuzungen als Knoten) für BAS-Klassifikation, Längen, Wurzeleffizienz und Kosten

### Vorprüfungen (Phase 0)
- Seed-Rauschen: 50 Seeds am Nominalpunkt (CV pro Output, Streuung zwischen Seeds als Ersatz für die Streuung zwischen Platten)
- Numerische Konvergenz: dx 0.01/0.005/0.0025 cm, dt 1/0.5/0.25 h (je 10 gepaarte Seeds)
- Strukturannahmen: Schichtdicke 0.02/0.005/0.1 cm, Mindestpfadlänge der Anastomose ×3/1/1.5/6
- Common Random Numbers: gepaarte Differenzen mit gleichen vs. unabhängigen Seeds
- 20 Zufallspunkte im Faktorraum: Abbruchrate und Laufzeit
- Langzeittest 180 h mit Wurzeln an vs. aus (Konvergenz bei 120 h): Zufluss neuer Hyphen, Einfluss auf die Outputs
- Ringauflösung 15 vs. 5 Ringe

### Sensitivitätsanalyse
- Elementareffekte nach Morris (1991), μ* nach Campolongo et al. (2007), Umsetzung mit SALib
- 22 Faktoren: alle kalibrierbaren Hyphenparameter, Variationskoeffizient der Spitzengeschwindigkeit, ein gemeinsamer Faktor
  auf alle übrigen Streuungen, die Dicke der Hyphenschicht, ein Kontrollfaktor ohne Wirkung
- Bereiche = Suchbereiche der Kalibrierung; log-gleichverteilt für Parameter über Größenordnungen
- 4 Stufen, 20 Trajektorien (aus 300 Kandidaten nach maximaler Streuung), 5 Replikate mit gemeinsamen Seeds entlang einer
  Trajektorie → 2300 Läufe
- Outputs:
  - pro Lauf, gemittelt über die Replikate: Längen, Spitzen, Frontposition und -geschwindigkeit, BAS-Anteil, Anastomosen,
    Wurzeleffizienz, Spitzengeschwindigkeiten, Winkel, Persistenzlänge, Zustand bei t*
  - gepoolt pro Design-Punkt: fitbasierte Ziele (alpha, beta, K1, K2, Frontgeschwindigkeit, Wellendauer)
- Wirksam: μ* über dem 95 %-Quantil der μ*-Verteilung bei reinem Seed-Rauschen (Monte Carlo aus den Standardfehlern)
- Zusätzlich: μ* skaliert auf die Breite des Literatur-Priors; Ausfallkarte; Korrelation der Ziele; Korrelation mit dem
  Zustand bei t*
- Separate Inokulum-Analyse (Wurzelphase, Wurzelarchitektur, Lage, Entstehungsdichte; 10 Trajektorien, 3 Replikate)

### Kalibrierung
- Gesucht ist x mit F(x) = 0, F_i = log(y_i/y*_i) / log(1+tol_i) (linear für Größen, die ≤ 0 werden können); |F_i| ≤ 1 heißt
  innerhalb der Toleranz.
- 21 Ziele aus Oyarte Galvez et al. 2025: Anastomose- und Verzweigungskonstante, Sättigungs- und Spitzendichte,
  Wellen- und Spitzengeschwindigkeit, Verzweigungswinkel, Persistenzlänge (Untergrenze) und die Zeitverläufe (20/40/60 h) von
  Runner-Länge, BAS-Anteil, Anastomosen und Wurzeleffizienz.
- Toleranzen aus der Streuung der 12 Replikate bzw. der Abbildungen; Zeitreihen um √3 erweitert (korrelierte Zeitpunkte);
  Ziele aus einer einzigen Kolonie um die Streuung zwischen Seeds erweitert.
- Bis zu 19 Parameter; Auswahl der freien Parameter nach der Morris-Analyse, übrige auf Nominalwert.
- Priors als zusätzliche Residuen: Literaturbereiche als ±2σ, sonst schwach um den Nominalwert (σ = 0.5 × Suchbereich); kein
  Prior auf Größen, deren Messung schon Ziel ist.
- Löser: scipy least_squares (Trust Region Reflective mit Grenzen) im normierten Parameterraum; Jacobi-Matrix aus
  Vorwärtsdifferenzen (Schritt 0.15), Differenzpunkte parallel.
- 16 Seeds pro Parameterpunkt mit gemeinsamen Seeds; Ziele über die Seeds gepoolt, Standardfehler per Bootstrap (200).
- Ergebnis: bester Punkt, Datenanteil pro Parameter, linearisierte ±1σ, Warnungen.
- Validierung mit nicht kalibrierten Größen (globale Effizienz, Kosten), danach varianzbasierte SA (Sobol über Ersatzmodell)
  für die kalibrierten Parameter.
