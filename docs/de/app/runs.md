---
title: Läufe in der Web-App und der Desktop-App
description: Wo die Läufe leben, wie lange sie bleiben und wie man sie behält, löscht und wiederherstellt.
sourceHash: 07a4256fedc9e26d
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

Der Browser behält zwei Dinge: Ihre Einstellungen der Oberfläche (etwa das Farbschema und die Breiten der Seitenleiste und des Optionsbereichs), die alle Tabs teilen, und das unfertige Formular jedes Tabs, das ein Neuladen dieses Tabs übersteht. Nach einem Neuladen behält das Formular seine Einstellungen und seine Beispieldateien, aber Ihre eigenen Dateien gehörten zur vorigen Seite: wählen Sie sie erneut.

Solange der Tab Läufe enthält, ein Lauf rechnet, ein Import eines Laufpakets offen ist oder das Formular Dateien enthält, die kein Lauf verwendet, fragt der Browser, bevor er den Tab schließt oder neu lädt. Die Liste zeigt die Größe jedes Laufs, weil alle Läufe eines Tabs den Speicher des Tabs teilen. Löschen Sie die Läufe, die Sie nicht mehr brauchen, um den Speicher freizugeben.

Wenn die Berechnung des Tabs fehlschlägt, zum Beispiel weil der Speicher ausgeht, sagt ein Banner oben im Hauptbereich, dass jeder Lauf des Tabs verloren ist und dass Laufpakete Läufe behalten. Die App startet die Berechnung selbst neu; wenn sie nicht starten kann, bietet das Banner an, die Seite neu zu laden. Wenn eine neue Version der Web-App erscheint, während der Tab solche Arbeit enthält, bietet ein Hinweis an, die Seite neu zu laden, statt sie selbst neu zu laden.

Ein Verweis auf einen Lauf, der nicht im Tab ist, etwa nach einem Neuladen, öffnet die Seite „This run is not in this tab“ mit ihrer Ursache und den Aktionen „Open the start page“ und „Import a run package“.

## Desktop-App

Die Desktop-App bewahrt jeden Lauf in einem Ordner ihres Laufordners auf, und die Liste zeigt jeden Lauf dieses Ordners, auch nach einem Neustart. „Change runs folder“ in der Befehlspalette oder die Ordner-Schaltfläche in der Kopfzeile wählt einen anderen Laufordner. Ein Laufpaket bringt einen Lauf auch auf einen anderen Computer.

Wenn die Berechnung der App fehlschlägt, sagt ein Banner, dass die rechnenden Läufe angehalten wurden und die beendeten Läufe im Laufordner bleiben; wenn die Berechnung nicht starten kann, bietet das Banner an, sie neu zu starten. Ein Verweis auf einen Lauf, der nicht im Laufordner ist, öffnet die Seite „This run is not in the runs folder“ mit denselben Aktionen wie in der Web-App.

## Läufe löschen

„Delete“ im Menü eines Laufs oder in der Fußzeile der Seitenleiste für ausgewählte Läufe entfernt die Läufe sofort aus der Liste und zeigt einen Hinweis mit „Undo“. „Undo“ bringt dieselben Läufe zurück. Wenn der Hinweis sich schließt oder wenn Sie die Seite verlassen oder verbergen, löscht die App die Läufe und gibt ihren Speicher oder ihre Ordner frei. Ein Lauf, der noch rechnet, fragt nach einer Bestätigung, und sein Löschen bricht die Berechnung zuerst ab.

Eine Variante funktioniert nach dem Löschen ihres Quelllaufs weiter, weil die Eingabedateien bleiben, solange ein Lauf sie verwendet. Wenn das offene Formular Dateien eines gelöschten Laufs verwendet, warnt die App vor dem Löschen, und das Formular verwirft diese Dateien.
