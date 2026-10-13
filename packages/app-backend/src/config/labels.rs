use itertools::Itertools;

const LABELS: [(&str, &str); 13] = [
  ("clock_std_dev", "Clock rate std. dev."),
  ("coalescent", "Coalescent time scale Tc"),
  ("coalescent_opt", "Optimize Tc"),
  ("coalescent_skyline", "Skyline coalescent"),
  ("confidence", "Date confidence intervals"),
  ("covariation", "Covariation-aware regression"),
  ("max_iter", "Iterations"),
  ("method_anc", "Ancestral method"),
  ("model", "Substitution model"),
  ("output_nwk_style", "Newick style"),
  ("output_selection", "Output files"),
  ("relax", "Relaxed clock (slack, coupling)"),
  ("reroot", "Reroot method"),
];

const WORDS: [(&str, &str); 22] = [
  ("aa", "amino-acid"),
  ("anc", "ancestral"),
  ("aln", "alignment"),
  ("alns", "alignments"),
  ("cdses", "CDSes"),
  ("csv", "CSV"),
  ("dot", "DOT"),
  ("dp", "DP"),
  ("gtr", "GTR"),
  ("iqd", "IQD"),
  ("json", "JSON"),
  ("mat", "MAT"),
  ("n", "number of"),
  ("nuc", "nucleotide"),
  ("nwk", "Newick"),
  ("opt", "optimization"),
  ("pb", "protobuf"),
  ("pc", "pseudocount"),
  ("pi", "pi"),
  ("rtt", "root-to-tip"),
  ("tsv", "TSV"),
  ("vcf", "VCF"),
];

pub fn setting_label(key_path: &[String]) -> String {
  key_path.iter().map(|part| key_label(part)).join(": ")
}

fn key_label(key: &str) -> String {
  if let Some((_, label)) = LABELS.iter().find(|(known, _)| *known == key) {
    return (*label).to_owned();
  }
  let words = key
    .split('_')
    .map(|word| {
      WORDS
        .iter()
        .find(|(short, _)| *short == word)
        .map_or(word, |(_, long)| long)
    })
    .join(" ");
  let mut chars = words.chars();
  chars
    .next()
    .map(|first| first.to_uppercase().chain(chars).collect())
    .unwrap_or_default()
}
