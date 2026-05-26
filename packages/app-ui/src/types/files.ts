export type FileSlotKind = "tree" | "alignment" | "dates" | "states" | "weights" | "vcfReference";

export interface LoadedFile {
  name: string;
  size: number;
}

export interface FileSlotConfig {
  kind: FileSlotKind;
  label: string;
  accept: string;
  description: string;
}

export const FILE_SLOTS: ReadonlyArray<FileSlotConfig> = [
  { kind: "tree", label: "Tree", accept: ".nwk,.nexus,.nex,.nhx", description: "Newick or Nexus format" },
  { kind: "alignment", label: "Alignment", accept: ".fasta,.fa,.fasta.xz,.fa.xz", description: "FASTA format" },
  { kind: "dates", label: "Dates", accept: ".csv,.tsv,.txt", description: "CSV/TSV with sampling dates" },
  { kind: "states", label: "States", accept: ".csv,.tsv,.txt", description: "CSV/TSV with discrete states" },
  { kind: "weights", label: "Weights", accept: ".csv,.tsv,.txt", description: "Equilibrium state weights" },
  { kind: "vcfReference", label: "VCF Reference", accept: ".fasta,.fa", description: "FASTA reference for VCF" },
];
