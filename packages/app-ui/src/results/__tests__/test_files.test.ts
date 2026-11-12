import type { OutputSelection, RunFile } from "@neherlab/app-contracts";

import { type DownloadEntry, downloadEntries } from "../files";

describe("download menu entries", () => {
  test.each([
    {
      analysis: "time tree",
      files: [
        file("timetree.auspice.json", "Auspice JSON", "auspice", true),
        file("timetree.nwk", "Dated tree (Newick, branch lengths in years)", "nwk", true),
        file("timetree.nexus", "Dated tree with node annotations (Nexus)", "nexus", true),
        file("timetree.rtt.tsv", "Root-to-tip table (TSV)", "rtt", false),
        file("timetree.trace.tsv", "Convergence trace (TSV)", "tracelog", false),
      ],
      expected: [
        "Run package",
        "Auspice JSON",
        "Dated tree (Newick, branch lengths in years)",
        "Dated tree with node annotations (Nexus)",
        "All files",
      ],
    },
    {
      analysis: "clock",
      files: [
        file("clock.rtt.tsv", "Root-to-tip table (TSV)", "rtt", true),
        file("clock.nwk", "Rerooted tree (Newick)", "nwk", true),
        file("clock.auspice.json", "Auspice JSON", "auspice", false),
      ],
      expected: ["Run package", "Root-to-tip table (TSV)", "Rerooted tree (Newick)", "All files"],
    },
    {
      analysis: "ancestral",
      files: [
        file("ancestral.fasta", "Ancestral sequences (FASTA)", "sequences", true),
        file("ancestral.nexus", "Annotated tree (Nexus)", "nexus", true),
        file("ancestral.auspice.json", "Auspice JSON", "auspice", true),
        file("ancestral.model.json", "Substitution model (JSON)", "model", false),
      ],
      expected: ["Run package", "Ancestral sequences (FASTA)", "Annotated tree (Nexus)", "Auspice JSON", "All files"],
    },
    {
      analysis: "a run without outputs",
      files: [],
      expected: ["Run package", "All files"],
    },
  ])("$analysis: the run package, the main formats in the given order, then all files", ({ files, expected }) => {
    expect(downloadEntries(files).map(entryText)).toStrictEqual(expected);
  });
});

function entryText(entry: DownloadEntry): string {
  if (entry.kind === "package") {
    return "Run package";
  }

  return entry.kind === "file" ? entry.file.label : "All files";
}

function file(path: string, label: string, kind: OutputSelection, main: boolean): RunFile {
  return { path, download_name: `run-1a2b3c4d.${path}`, label, size: 1, kind, description: "", main };
}
