---
title: Arbeitsbereich
description: Der Optionsbereich, die Fakten der Eingaben, die Einstellungen, die Prüfungen und die Vorschau vor dem ersten Lauf.
sourceHash: df681bb662ed270d
---

Der Arbeitsbereich zeigt links die Optionen einer Analyse und rechts den Hauptbereich. Dasselbe Layout dient einer neuen Analyse und einem Lauf: während Sie eine Annahme ändern, bleibt das Ergebnis, das zu der Änderung führte, sichtbar. Ziehen Sie die Grenze zwischen den beiden Bereichen, um ihre Breiten zu ändern. Auf einem schmalen Bildschirm öffnen sich die Optionen in einem Einschub, und die Schaltfläche Run bleibt in einer Leiste am unteren Rand.

Die Kopfzeile des Arbeitsbereichs nennt den Datensatz: den Namen des Beispiels oder die Dateinamen Ihrer Daten.

## Optionsbereich

Der Optionsbereich hat vier Teile, von oben nach unten.

### Analyse

Wählen Sie die Analyse: Time tree (standardmäßig ausgewählt), Clock signal, Ancestral sequences, Optimize tree oder Prune tree. Die Info-Schaltfläche neben der Auswahl beschreibt die Analyse und nennt ihren `treetime`-Befehl.

### Daten

Die Eingabefelder (Baum, Alignment, Metadaten) mit ihren Dateien und die Fakten, die TreeTime in ihnen erkannt hat: die Zahl der Blätter, die Länge und das Alphabet des Alignments, das Trennzeichen der Metadaten, die ID-Spalte und die Datumsspalte mit der Zahl der gelesenen Daten. Zum Beispiel: „Dates: column `date`, 100 of 100 read as YYYY-MM-DD“. „Change“ neben einem Fakt öffnet die Einstellung, die die Erkennung ersetzt.

### Haupteinstellungen

Die wichtigsten Annahmen der Analyse, etwa die Uhrrate, die Wurzel, der Uhrfilter, der Koaleszenz-Prior und die Konfidenzintervalle. Jede Zeile zeigt den Wert, den der Lauf verwenden wird, auch wenn Sie nichts setzen, etwa „Automatic: least-squares rerooting“. Die Hilfe-Schaltfläche einer Zeile erklärt die Einstellung und ihre Optionen. Der Abschnittstitel zeigt, wie viele Haupteinstellungen von ihren Standardwerten abweichen.

### Alle Einstellungen

„All settings“ öffnet jede Einstellung der Analyse in einem Bereich rechts, nach Zweck gruppiert:

- **Input data**: wie TreeTime die Eingaben liest
- **Analysis**: das Modell und die Inferenz, nach Thema gruppiert
- **Output format**: die Einstellungen, die den Inhalt der Ausgabedateien ändern, etwa die Annotationen der Newick-Ausgabe
- **Tree ordering and reproducibility**: die Reihenfolge der Knoten im Ausgabebaum und der Startwert (Seed) der Zufallsschritte

Die Suche findet eine Einstellung nach ihrem Namen, ihrer Beschreibung oder ihrer Befehlszeilenoption. „Changed only“ zeigt nur die Einstellungen, die von ihren Standardwerten abweichen, und die Rücksetz-Schaltfläche einer Einstellung setzt sie auf ihren Standardwert zurück. Die Hilfe jeder Einstellung nennt ihre Befehlszeilenoption, etwa „Command line: `--clock-filter 4`“. Die App setzt die Ausgabepfade jedes Laufs selbst, daher stehen sie nicht im Formular.

### Prüfungen und Run

Die Fußzeile unten im Bereich prüft die Optionen, während Sie sie bearbeiten. Sie zeigt „Ready to run.“ oder die Befunde: einen Befund, der den Lauf blockiert, eine Warnung oder einen Rat. Ein Befund nennt die Einstellungen mit ihren Bezeichnungen und kann eine Korrektur anbieten: eine Schaltfläche, die die Einstellungen ändert, etwa „Use covariation“. Run startet die Analyse; das Tastenkürzel ist Strg+Eingabe (Cmd+Eingabe unter macOS).

## Optionsmenü

Die Menü-Schaltfläche des Optionsbereichs bietet:

- **Command line**: die Optionen als `treetime`-Befehl und als YAML-Konfiguration; siehe [Befehlszeile](/app/export#befehlszeile). Das Tastenkürzel ist C
- **Apply settings from a run**: die Einstellungen, die ein gewählter Lauf gegenüber den Standardwerten geändert hat, ersetzen die Einstellungen der Optionen. Die Eingaben bleiben. Die Liste zeigt die beendeten Läufe derselben Analyse
- **Load config (YAML)**: die Einstellungen einer YAML-Konfigurationsdatei laden, etwa einer, die `treetime` mit `--config` liest

## Vor dem ersten Lauf

Vor dem Lauf zeigt der Hauptbereich, was TreeTime aus den Eingaben ablesen kann:

- **Input data**: die Fakten des Baums, des Alignments und der Metadaten
- **Tree**: der Eingabebaum ohne Daten, in Divergenzeinheiten
- **Root-to-tip regression**: die Vorschau der Uhr, der Abstand jeder datierten Probe von der Wurzel gegen ihr Datum, mit der Regressionsgeraden und den Ausreißern des Uhrfilters in Rot. Sie verwendet die Einstellungen der Wurzel und des Uhrfilters aus den Optionen, sodass Sie das zeitliche Signal vor dem Lauf prüfen können

Die Vorschau berechnet Bäume mit bis zu 5100 Blättern sofort. Für einen größeren Baum berechnet „Update preview“ sie auf Anfrage.

## Einstellungen eines gezeigten Laufs

Ein Klick auf einen Lauf in der Liste zeigt seine Ergebnisse und seine Optionen. Die Optionen sind schreibgeschützt, bis Sie eine ändern: die erste Änderung macht sie zu einer [Variante](/app/variants) des Laufs. Eine Leiste über den Ergebnissen listet dann die Einstellungen, die vom gezeigten Lauf abweichen, etwa „Settings changed since this run: coalescent prior, confidence. Run to update.“ Run startet die Variante, und „Revert to this run's settings“ stellt die Optionen des Laufs wieder her.

„New analysis“ in der Seitenleiste beginnt mit den Standardeinstellungen und ohne Daten. „Reset form“ leert das Formular, und „Undo“ in seinem Hinweis stellt das Formular wieder her.
