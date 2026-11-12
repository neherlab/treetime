---
title: Läufe in der Web-App und der Desktop-App
description: Wo die Läufe leben, wie lange sie bleiben und wie man sie behält, löscht und wiederherstellt.
sourceHash: d440f5ca00f99eb5
---

Die Web-App und die Desktop-App bewahren ihre Läufe an verschiedenen Orten auf, daher bestehen ihre Lauflisten verschieden lange.

|                                             | Web-App                                               | Desktop-App                       |
| ------------------------------------------- | ----------------------------------------------------- | --------------------------------- |
| Läufe                                       | Im Speicher des Browser-Tabs                          | Im Laufordner auf der Festplatte  |
| Liste                                       | Die Läufe dieses Tabs                                 | Jeder Lauf des Laufordners        |
| Nach einem Neuladen oder Schließen des Tabs | Die Läufe sind weg                                    | Die Läufe bleiben                 |
| Einen Lauf behalten                         | Sein [Laufpaket](/app/export#laufpaket) herunterladen | Der Laufordner behält ihn         |
| Unfertiges Formular                         | Für diesen Tab über ein Neuladen hinweg behalten      | In der Einstellungsdatei behalten |

## Web-App

Die Web-App rechnet im Browser-Tab und speichert keinen Lauf im Browser. Ein Lauf besteht so lange wie der Tab: ein Neuladen oder ein geschlossener Tab entfernt jeden Lauf des Tabs. Um einen Lauf zu behalten oder ihn in einem anderen Tab, einem anderen Browser oder der Desktop-App zu öffnen, laden Sie sein Laufpaket herunter und importieren es dort. Jeder Tab hat seine eigenen Läufe, und zwei Tabs der Web-App sehen die Läufe des anderen nicht.

Der Browser behält zwei Dinge: Ihre Einstellungen der Oberfläche (etwa das Farbschema und die Breite der Seitenleiste), die alle Tabs teilen, und das unfertige Formular jedes Tabs, das ein Neuladen dieses Tabs übersteht. Nach einem Neuladen behält das Formular seine Einstellungen und seine Beispieldateien, aber Ihre eigenen Dateien gehörten zur vorigen Seite: wählen Sie sie erneut.

Solange der Tab Läufe enthält, fragt der Browser, bevor er den Tab schließt oder neu lädt. Die Liste zeigt die Größe jedes Laufs, weil alle Läufe eines Tabs den Speicher des Tabs teilen. Löschen Sie die Läufe, die Sie nicht mehr brauchen, um den Speicher freizugeben.

Wenn die Berechnung des Tabs fehlschlägt, zum Beispiel weil der Speicher ausgeht, sagt ein Banner oben im Hauptbereich, dass jeder Lauf des Tabs verloren ist und dass Laufpakete Läufe behalten. Die App startet die Berechnung selbst neu; wenn sie nicht starten kann, bietet das Banner an, die Seite neu zu laden. Wenn eine neue Version der Web-App erscheint, während ein Tab Läufe enthält, bietet ein Hinweis an, die Seite neu zu laden, statt sie selbst neu zu laden.

Ein Verweis auf einen Lauf, der nicht im Tab ist, etwa nach einem Neuladen, öffnet eine Seite, die das sagt, mit Aktionen, die die Laufliste öffnen oder ein Laufpaket importieren.

## Desktop-App

Die Desktop-App bewahrt jeden Lauf in einem Ordner ihres Laufordners auf, und die Liste zeigt jeden Lauf dieses Ordners, auch nach einem Neustart. „Change runs folder“ in der Befehlspalette oder die Ordner-Schaltfläche in der Kopfzeile wählt einen anderen Laufordner. Ein Laufpaket bringt einen Lauf auch auf einen anderen Computer.

## Läufe löschen

„Delete“ im Menü eines Laufs oder in der Fußzeile der Seitenleiste für ausgewählte Läufe entfernt die Läufe sofort aus der Liste und zeigt einen Hinweis mit „Undo“. „Undo“ bringt dieselben Läufe zurück. Wenn der Hinweis sich schließt, löscht die App die Läufe und gibt ihren Speicher oder ihre Ordner frei. Ein Lauf, der noch rechnet, fragt nach einer Bestätigung, und sein Löschen bricht die Berechnung zuerst ab.

Eine Variante funktioniert nach dem Löschen ihres Quelllaufs weiter, weil die Eingabedateien bleiben, solange ein Lauf sie verwendet. Wenn das offene Formular Dateien eines gelöschten Laufs verwendet, warnt die App vor dem Löschen, und das Formular verwirft diese Dateien.
