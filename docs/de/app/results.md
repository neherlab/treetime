---
title: Ergebnisse
description: Die Kopfzeile des Laufs, die Übersicht und die Ergebnisansichten jeder Analyse.
sourceHash: 174e1b9668e8e4c9
---

Eine Laufseite zeigt das Ergebnis eines Laufs. Die Optionen des Laufs bleiben im Bereich links, sodass Sie eine Annahme neben dem Ergebnis ändern können; siehe [Einstellungen eines gezeigten Laufs](/app/workspace#einstellungen-eines-gezeigten-laufs).

## Kopfzeile des Laufs

Die Kopfzeile zeigt den Namen des Laufs, den Status, die Analyse, den Datensatz (den Namen des Beispiels, sonst den Dateinamen des Eingabebaums), das Datum des Laufs und seine Dauer. Die Stift-Schaltfläche neben dem Namen benennt den Lauf um. Die Aktionen der Kopfzeile:

- **Download**: das Laufpaket und die wichtigsten Ausgabedateien; siehe [Download-Menü](/app/export#download-menü)
- **Command line**: der Befehl und die Konfiguration, die den Lauf wiederholen; siehe [Befehlszeile](/app/export#befehlszeile)
- **Create variant**: ein neuer Lauf aus den Einstellungen dieses Laufs; siehe [Varianten](/app/variants)
- **Menü**: umbenennen, anheften, mit der Quelle vergleichen, Laufpaket herunterladen und löschen

Wenn der Lauf Warnungen geschrieben hat, nennt eine Zeile unter der Kopfzeile ihre Zahl und die erste Ursache, etwa „3 warnings: 3 tips dated outside their input date range“. „Show“ listet sie alle. Eine Warnung steht auch neben der Schätzung, die sie betrifft.

## Ansichten

Die Reiter unter der Kopfzeile wählen die Ansicht. Jede Ansicht exportiert ihre eigenen Daten. Die Auswahl eines Knotens und die Filter des Baums bleiben, wenn Sie die Ansicht wechseln.

| Ansicht    | Inhalt                                                                                                                                           | Analysen                                      |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------- |
| Overview   | Die wichtigsten Schätzungen und der Baum                                                                                                         | Alle                                          |
| Clock      | Die Wurzel-Blatt-Regression, die Residuen, die Regressionsstatistiken, die Uhr-Ausreißer und die undatierten Proben                              | Time tree                                     |
| Population | Der Koaleszenz-Prior: die Skyline oder die Koaleszenz-Zeitskala, und die Koaleszenz-Einstellungen                                                | Time tree                                     |
| Model      | Das Substitutionsmodell und die Entropie jeder Stelle                                                                                            | Time tree, Ancestral sequences, Optimize tree |
| Log        | Das Protokoll des Laufs, nach Abschnitten, mit Filter, Suche und Kopieraktion                                                                    | Alle                                          |
| Settings   | Die Einstellungen des Laufs, mit einem Download seiner YAML-Konfiguration; „Changed only“ zeigt nur die Einstellungen, die der Lauf geändert hat | Alle                                          |

Ein Lauf, der nicht beendet ist, zeigt nur Overview, Log und Settings.

### Übersicht eines Zeitbaums

Der Schätzungsbereich zeigt das Datum der Wurzel mit seinem Intervall und die Uhrrate mit ihrer Standardabweichung, etwa „7.6e-4 per site per year, ± 4.1e-5“, dann das zeitliche Signal und die Zahl der verwendeten Proben. „Run details“ zeigt das Uhrmodell, den Koaleszenz-Prior, die Log-Likelihood, die Zahl der Iterationen und die Dauer.

Wenn der Lauf kein Intervall für die Daten berechnet hat, sagt der Bereich „No interval computed“. Wenn der Lauf die Konfidenzintervalle ausgeschaltet hatte, bietet er außerdem „Compute the interval“ an: das startet eine Variante des Laufs mit eingeschalteten Konfidenzintervallen.

Der datierte Baum ist das größte Element. Auf einem breiten Bildschirm steht das Wurzel-Blatt-Diagramm daneben. Ein Klick auf einen Knoten im Baum oder auf einen Punkt im Diagramm wählt die Probe in beiden aus, und ein Seitenbereich zeigt ihre Details: ihr Datum, ihre Uhrdaten, ihre Mutationen und, falls die Uhranpassung sie ausgeschlossen hat, den Grund.

### Übersicht der anderen Analysen

- **Clock signal**: das Wurzel-Blatt-Diagramm mit dem Baum daneben und die Schätzungen der Regression
- **Ancestral sequences**: der Baum, gefärbt nach den Mutationen an der Stelle, die auf den meisten Ästen mutiert, die Mutationen jedes Asts und die Stellen, die auf mehreren Ästen mutiert sind
- **Optimize tree** und **Prune tree**: der Ergebnisbaum und eine Tabelle mit Statistiken des Eingabebaums und des Ergebnisbaums

### Ansicht Clock

Die Ansicht Clock zeigt die Regressionsstatistiken so, wie Tempest sie listet: die Uhrrate, das Datum der Wurzel aus dem x-Achsenabschnitt, den Korrelationskoeffizienten r, R² und das mittlere Residuenquadrat. Die Info-Schaltfläche jeder Statistik gibt ihre Definition. Die Listen der Uhr-Ausreißer und der undatierten Proben verweisen auf den Baum und das Diagramm.

## Baumansicht

Eine Werkzeugleiste über dem Baum bietet die Färbung, das Layout, die Suche und den Zoom auf die Wurzel. Die Suche listet die Proben, deren Namen zu Ihrer Eingabe passen, auch bei kleinen Tippfehlern, und wählt die Probe aus, die Sie anklicken. „View options“ öffnet alle Optionen der Auspice-Ansicht, darunter den Download des Baums, wie er gezeigt wird.

## Laufende und fehlgeschlagene Läufe

Während ein Lauf rechnet, zeigt der Hauptbereich seine Abschnitte, die Uhrrate jeder Iteration und eine Abbruchaktion. Die Optionen bleiben bearbeitbar, sodass Sie den nächsten Lauf vorbereiten können.

Ein fehlgeschlagener Lauf zeigt den Fehler mit seinen Ursachen und den betroffenen Proben sowie die Abschnitte des Protokolls. „Edit and run again“ öffnet die Optionen an der Einstellung, die der Fehler nennt.
