const WORDS = {
  aa: "amino-acid",
  anc: "ancestral",
  aln: "alignment",
  alns: "alignments",
  cdses: "CDSes",
  csv: "CSV",
  dot: "DOT",
  dp: "DP",
  gtr: "GTR",
  iqd: "IQD",
  json: "JSON",
  mat: "MAT",
  n: "number of",
  nuc: "nucleotide",
  nwk: "Newick",
  opt: "optimization",
  pb: "protobuf",
  pc: "pseudocount",
  pi: "pi",
  rtt: "root-to-tip",
  tsv: "TSV",
  vcf: "VCF",
} satisfies Record<string, string>;

const LABELS = {
  clock_std_dev: "Clock rate std. dev.",
  coalescent: "Coalescent time scale Tc",
  coalescent_opt: "Optimize Tc",
  coalescent_skyline: "Skyline coalescent",
  confidence: "Date confidence intervals",
  covariation: "Covariation-aware regression",
  max_iter: "Iterations",
  method_anc: "Ancestral method",
  model: "Substitution model",
  output_nwk_style: "Newick style",
  output_selection: "Output files",
  relax: "Relaxed clock (slack, coupling)",
  reroot: "Reroot method",
} satisfies Record<string, string>;

const WORD_MAP = new Map<string, string>(Object.entries(WORDS));

const LABEL_MAP = new Map<string, string>(Object.entries(LABELS));

export function settingLabel(key: string): string {
  const path = key.split(".");

  if (path.length > 1) {
    return path.map((part) => settingLabel(part)).join(": ");
  }

  const known = LABEL_MAP.get(key);

  if (known !== undefined) {
    return known;
  }

  const words = key.split("_").map((word) => WORD_MAP.get(word) ?? word);
  const [first = "", ...rest] = words;

  return [first.charAt(0).toUpperCase() + first.slice(1), ...rest].join(" ");
}
