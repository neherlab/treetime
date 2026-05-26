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

export const EXAMPLE_DATASETS: ReadonlyArray<string> = [
  "dengue/100",
  "dengue/1000",
  "dengue/20",
  "dengue/2000",
  "dengue/500",
  "ebola/100",
  "ebola/20",
  "ebola/362",
  "flu/h3n2/20",
  "flu/h3n2/200",
  "flu/h3n2/500",
  "lassa/L/20",
  "lassa/L/200",
  "lassa/L/50",
  "lassa/L/500",
  "mpox/clade-ii/100",
  "mpox/clade-ii/1000",
  "mpox/clade-ii/20",
  "mpox/clade-ii/2000",
  "mpox/clade-ii/500",
  "rsv/a/100",
  "rsv/a/1000",
  "rsv/a/20",
  "rsv/a/2000",
  "rsv/a/500",
  "sc2/2844",
  "tb/100",
  "tb/149",
  "tb/20",
  "zika/20",
  "zika/86",
];
