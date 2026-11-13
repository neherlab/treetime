---
title: Die TreeTime-App
description: Einen Datensatz öffnen, sein Uhrsignal prüfen, eine Analyse ausführen und das Ergebnis in der Web-App oder der Desktop-App lesen.
sourceHash: 3b4a6c6d05f978a5
---

Die TreeTime-App führt die Analysen des Befehlszeilenprogramms `treetime` mit Ihren eigenen Daten oder mit einem Beispieldatensatz aus. Zwei Apps teilen eine Oberfläche:

- **Web-App** ([treetime.neherlab.click](https://treetime.neherlab.click)): TreeTime läuft im Browser-Tab. Ihre Dateien bleiben im Browser, und nichts wird auf einen Server hochgeladen. Die Läufe leben im Tab; siehe [Läufe in der Web-App und der Desktop-App](/app/runs)
- **Desktop-App**: TreeTime läuft auf Ihrem Computer und behält jeden Lauf in einem Laufordner auf der Festplatte

Beide Apps rechnen mit demselben Code wie das Befehlszeilenprogramm, und jeder Lauf gibt seine [Befehlszeile](/app/export#befehlszeile) an, sodass Sie ihn mit `treetime` wiederholen können.

Der übliche Arbeitsablauf:

1. Einen Datensatz öffnen: Ihre Dateien ablegen oder ein Beispiel öffnen
2. Die Fakten der Eingaben und die Vorschau der Uhr ansehen und korrigieren, was die App erkannt hat
3. Die Analyse ausführen und die Übersicht lesen
4. Eine Annahme in einer [Variante](/app/variants) ändern und den Unterschied zum Quelllauf lesen
5. Das Ergebnis [exportieren](/app/export) oder das Laufpaket herunterladen, um den Lauf zu behalten

## Startseite

Die Startseite öffnet sich, wenn Sie die App unter ihrer Adresse öffnen.

- **Ablagefläche**: legen Sie hier den Baum, das Alignment und die Metadaten ab, oder klicken Sie auf „Choose files“. Die App liest den Anfang jeder Datei, bevor etwas hochgeladen wird: eine Zip-Datei wird als Laufpaket importiert, eine Datei, die sich als TreeTime-Konfiguration lesen lässt, wird in die Optionen geladen, und jede andere Datei ist eine Eingabe. Sie liest den Inhalt jeder Eingabe, um ihre Art zu bestimmen: die Dateiendung ist nur ein Hinweis. Sie legt jede Eingabe in ihr Feld und öffnet die Optionen einer neuen Analyse. Eine Datei, die die App nicht zuordnen kann, bleibt draußen, und ein Hinweis nennt für jede solche Datei den Grund. Unter der Schaltfläche listet die App die Formate, die jede Eingabe annimmt
- **Konfiguration**: eine YAML-Konfiguration auf der Ablagefläche lädt ihre Einstellungen in die Optionen, wie „Load config (YAML)“ im [Optionsmenü](/app/workspace#optionsmenü)
- **Laufpaket**: ein Laufpaket (`.zip`) auf der Ablagefläche importiert den Lauf mit seinen Ergebnissen; siehe [Laufpaket](/app/export#laufpaket)
- **Beispiele**: eine Karte pro Beispieldatensatz, mit dem Organismus und dem Segment, der Zahl der Proben, dem Datumsbereich und der Sequenzlänge. „Open“ öffnet den Datensatz in einer neuen Analyse, mit dem Zeitbaum (Time tree) ausgewählt. In der Desktop-App ohne Beispiele bietet die Galerie an, sie herunterzuladen
- **Letzte Arbeit**: die letzten Arbeitslinien, jede mit dem Namen, dem Status und der wichtigsten Schätzung ihres neuesten Laufs

Sie können Dateien auch an jeder anderen Stelle der App ablegen: sie kommen in die Eingabefelder der offenen Analyse.

## Laufliste

Die Seitenleiste listet die Läufe. Ein Lauf und die Varianten, die aus ihm entstanden, bilden eine Arbeitslinie: jede Variante steht unter ihrem Quelllauf. Ein Eintrag zeigt den Namen des Laufs, ein Statussymbol, die wichtigste Schätzung (etwa das Datum der Wurzel), den Tag des Laufs und seine Größe. „Filter by name“ und der Analysefilter schränken die Liste ein.

Das Menü eines Eintrags hat diese Aktionen:

- **Rename**: dem Lauf einen eigenen Namen geben. Die App benennt jeden Lauf nach dem Datensatz, der Analyse und den abweichenden Einstellungen, etwa „ebola/100, time tree: skyline, confidence“
- **Pin**: die Arbeitslinie des Laufs oben in der Liste halten
- **Compare with source**: den [Vergleich](/app/variants#zwei-läufe-vergleichen) mit dem Lauf öffnen, aus dem die Variante entstand
- **Download run package**: siehe [Laufpaket](/app/export#laufpaket)
- **Delete**: siehe [Läufe löschen](/app/runs#läufe-löschen)

Um mit mehreren Läufen zu arbeiten, wählen Sie „Select“ im Menü eines Eintrags, oder klicken Sie die Einträge mit Strg (Cmd unter macOS) oder Umschalt an. Die Fußzeile der Seitenleiste bietet dann „Compare“ für genau zwei beendete Läufe, mit dem zuerst ausgewählten Lauf als Lauf A, „Delete“ und „Clear“.

„New analysis“ oben in der Seitenleiste beginnt eine Analyse mit den Standardeinstellungen und ohne Daten. Die Befehlspalette (Strg+K, Cmd+K unter macOS) findet Läufe, Beispiele und Aktionen nach ihrem Namen.
