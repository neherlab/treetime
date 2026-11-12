---
title: Export
description: Ausgabedateien und Laufpakete herunterladen, ein Laufpaket importieren und einen Lauf mit dem Befehlszeilenprogramm wiederholen.
sourceHash: 0d28d205941c33b3
---

## Download-Menü

„Download“ in der Kopfzeile des Laufs listet, in dieser Reihenfolge:

- **Run package (.zip)**: die Eingaben, die Einstellungen, das Protokoll und die Ausgaben des Laufs; siehe [Laufpaket](#laufpaket)
- **Wichtigste Formate** der Analyse:
  - Time tree: Auspice-JSON, der datierte Baum in Newick (Astlängen in Jahren) und Nexus
  - Clock signal: die Wurzel-Blatt-Tabelle und der neu gewurzelte Baum in Newick
  - Ancestral sequences: die Vorfahrensequenzen in FASTA, der annotierte Baum in Nexus und Auspice-JSON
  - Optimize tree und Prune tree: der Baum in Newick
- **All files...**: ein Dialog mit jeder Ausgabedatei des Laufs und einer Beschreibung jeder Datei

Der Name einer heruntergeladenen Datei beginnt mit dem Namen des Laufs und den ersten Zeichen der Lauf-ID, sodass die Dateien verschiedener Läufe einander im Download-Ordner nicht ersetzen.

Wenn der Browser in eine gewählte Datei speichern kann (die Seite kommt über HTTPS oder von `localhost`, und der Browser hat den Speicherdialog der File System Access API, wie Chrome und Edge), fragt die Web-App, wohin sie speichern soll, und schreibt die Datei auf die Festplatte, während sie sie liest. Sonst hält die Web-App die ganze Datei im Speicher der Seite, bis der Download beginnt, daher kann der Download eines großen Laufpakets dort fehlschlagen.

## Laufpaket

Ein Laufpaket ist eine Zip-Datei, die einen ganzen Lauf enthält: seinen Datensatz mit den Einstellungen und der Zusammenfassung der Ergebnisse, das Protokoll, die Eingabedateien, die Ausgabedateien, die YAML-Konfiguration und ein Shell-Skript mit der Befehlszeile. Laden Sie es im Download-Menü oder im Menü des Laufs in der Liste herunter.

Ein Laufpaket ist der Weg, einen Lauf der Web-App zu behalten und einen Lauf zu teilen: legen Sie es auf der Startseite oder an einer beliebigen Stelle der App ab, um es zu importieren. Der importierte Lauf zeigt dieselben Ergebnisse, und Sie können Varianten davon erstellen, als hätten Sie ihn in dieser App ausgeführt. Die Desktop-App importiert die Pakete der Web-App und die Web-App die der Desktop-App. Die App prüft ein Paket vor dem Import: der Datensatz muss lesbar sein, der Lauf muss beendet sein, und jede Eingabedatei muss die Größe und die SHA-256-Prüfsumme haben, die der Datensatz angibt.

Das Paket wiederholt den Lauf auch mit dem Befehlszeilenprogramm. Entpacken Sie es und führen Sie das Skript im entpackten Ordner `treetime-<run-id>` aus:

```sh
unzip <package>.zip
cd treetime-<run-id>
sh command.sh
```

Das Skript führt `treetime` mit der Konfiguration `run.yaml` aus, die die Eingaben aus `inputs/` liest und die Ausgaben nach `out/` schreibt, mit dem Startwert des Laufs.

## Befehlszeile

„Command line“ (Tastenkürzel C) zeigt den Lauf oder die offenen Optionen als `treetime`-Befehl oder als YAML-Konfiguration, mit einer Kopier-Schaltfläche. Der Befehl enthält die Eingaben, die Einstellungen, die von den Standardwerten abweichen, die Ausgaben, die die App jedem Lauf hinzufügt, und den Startwert des Laufs, mit den Pfaden des Laufpakets. Der Dialog gibt außerdem:

- **Methods**: einen Absatz für den Methodenteil einer Veröffentlichung, mit der Analyse, der TreeTime-Version und den wichtigsten Einstellungen
- **Citation**: die Veröffentlichung, die für die Analyse zu zitieren ist
