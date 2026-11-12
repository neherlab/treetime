---
title: Varianten und Vergleich
description: Eine Annahme eines Laufs ändern, den Unterschied zum Quelllauf lesen und zwei Läufe vergleichen.
sourceHash: 0d4585aa2024da96
---

Eine Variante ist ein Lauf aus den Einstellungen eines anderen Laufs, ihrer Quelle, mit einer oder mehreren geänderten Einstellungen. Eine Variante und ihre Quelle bilden eine Arbeitslinie in der [Laufliste](/app#laufliste).

## Eine Variante erstellen

„Create variant“ in der Kopfzeile des Laufs macht die Optionen des Laufs bearbeitbar. Der Bereich nennt dann den Quelllauf und listet die Einstellungen, die Sie geändert haben, und die Liste wächst, während Sie bearbeiten. Sie können eine Variante auch mit der Änderung einer beliebigen Einstellung eines gezeigten Laufs beginnen; siehe [Einstellungen eines gezeigten Laufs](/app/workspace#einstellungen-eines-gezeigten-laufs). Run startet die Variante. Ihr Name listet die geänderten Einstellungen, etwa „ebola/100, time tree: skyline“.

Eine Variante behält den Startwert (Seed) ihres Quelllaufs. TreeTime löst Polytomien mit einem Zufallsprozess auf, daher können zwei Läufe mit verschiedenen Startwerten leicht verschiedene Schätzungen geben. Mit demselben Startwert kommt ein Unterschied zwischen der Variante und ihrer Quelle nur von den geänderten Einstellungen. Um einen neuen Startwert zu ziehen, leeren Sie den Startwert unter All settings, im Abschnitt „Tree ordering and reproducibility“.

## Unterschied zur Quelle

Die Übersicht einer Variante zeigt eine Zeile mit der Änderung der wichtigsten Schätzungen gegenüber dem Quelllauf, etwa „Root date -46 d, clock rate -7.9 % against the source“, mit einem Verweis auf den vollständigen Vergleich.

## Zwei Läufe vergleichen

Öffnen Sie den Vergleich mit „Compare with source“ im Menü einer Variante, oder wählen Sie zwei beendete Läufe in der Liste aus und klicken Sie auf „Compare“: der Lauf, den Sie zuerst ausgewählt haben, ist Lauf A. Die Tausch-Schaltfläche vertauscht Lauf A und Lauf B.

- **Estimates**: das Datum der Wurzel, die Uhrrate, die Breite des Intervalls, die Log-Likelihood, der Korrelationskoeffizient r und die ausgeschlossenen Proben beider Läufe, mit dem Unterschied in der Spalte „Difference (B - A)“. Der Vergleich markiert die Log-Likelihoods als nicht vergleichbar, wenn die Läufe sich in ihren Daten, ihrem Koaleszenz-Prior, ihrem Substitutionsmodell oder der Topologie ihrer Ergebnisbäume unterscheiden, weil die beiden Werte dann zu verschiedenen Modellen gehören
- **Trees**: ob die Bäume dieselbe Topologie haben und wie weit sich jeder gemeinsame Vorfahre zwischen den beiden Läufen in der Zeit verschiebt
- **Settings that differ**: die Einstellungen mit verschiedenen Werten in den beiden Läufen. Die Ausgabeeinstellungen, die die App für jeden Lauf setzt, stehen nicht in der Liste

Wenn die beiden Läufe verschiedene Analysen sind, sagt die Seite das.
