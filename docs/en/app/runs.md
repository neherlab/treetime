---
title: Runs in the web app and the desktop app
description: Where the runs live, how long they stay, and how to keep, delete, and recover them.
---

The web app and the desktop app keep their runs in different places, so their run histories last for different times.

|                                       | Web app                                             | Desktop app                  |
| ------------------------------------- | --------------------------------------------------- | ---------------------------- |
| Runs                                  | In the memory of the browser tab                    | In the runs folder on disk   |
| History                               | The runs of this tab                                | Every run of the runs folder |
| After a reload or when the tab closes | The runs are gone                                   | The runs stay                |
| Keep a run                            | Download its [run package](/app/export#run-package) | The runs folder keeps it     |
| Unfinished form                       | Kept for this tab through a reload                  | Kept in the settings file    |

## Web app

The web app computes in the browser tab and stores no run in the browser. A run stays as long as the tab: a reload or a closed tab removes every run of the tab. To keep a run, or to open it in another tab, another browser, or the desktop app, download its run package and import it there. Each tab has its own runs, and two tabs of the web app do not see the runs of each other.

The browser keeps two things: your preferences (such as the theme and the width of the sidebar), which all tabs share, and the unfinished form of each tab, which survives a reload of that tab. After a reload, the form keeps its settings and its example files, but your own files belonged to the previous page: choose them again.

While the tab holds runs, a run computes, a run package import is open, or the form holds files that no run uses, the browser asks before it closes or reloads the tab. The history shows the size of each run, because all runs of a tab share the memory of the tab. Delete the runs that you no longer need to free the memory.

When the computation of the tab fails, for example when it runs out of memory, a banner at the top of the main area says that every run of the tab is lost and that run packages keep runs. The app restarts the computation by itself; when it cannot start, the banner offers to reload the page. When a new version of the web app is published while the tab holds such work, a notice offers to reload the page instead of reloading it by itself.

A link to a run that is not in the tab, for example after a reload, opens the page "This run is not in this tab" with its cause and the actions "Open the run list" and "Import a run package".

## Desktop app

The desktop app keeps each run in a folder of its runs folder, and the history lists every run of that folder, also after a restart. "Change runs folder" in the command palette, or the folder button in the header, selects another runs folder. A run package also moves a run between computers.

When the computation of the app fails, a banner says that the runs that were computing stopped and that the finished runs stay in the runs folder; when the computation cannot start, the banner offers to restart it. A link to a run that is not in the runs folder opens the page "This run is not in the runs folder" with the same actions as in the web app.

## Delete runs

"Delete" in the menu of a run, or in the sidebar footer for selected runs, removes the runs from the history at once and shows a notice with Undo. Undo puts back the same runs. When the notice closes, the app deletes the runs and frees their memory or their folders. A run that still computes asks for confirmation, and its deletion cancels the computation first.

A variant keeps working after the deletion of its source run, because the input files stay as long as a run uses them. When the open form uses files of a deleted run, the app warns before the deletion, and the form drops those files.
