use itertools::Itertools;

const LABELS: [(&str, &str); 11] = [
  ("clock_std_dev", "Clock rate std. dev."),
  ("coalescent", "Coalescent prior"),
  ("confidence", "Date confidence intervals"),
  ("max_iter", "Iterations"),
  ("method_anc", "Ancestral method"),
  ("model", "Substitution model"),
  ("n_skyline", "Skyline grid points"),
  ("output_nwk_style", "Newick style"),
  ("output_selection", "Output files"),
  ("relax", "Relaxed clock (slack, coupling)"),
  ("reroot", "Reroot method"),
];

const WORDS: [(&str, &str); 19] = [
  ("aa", "amino-acid"),
  ("anc", "ancestral"),
  ("aln", "alignment"),
  ("alns", "alignments"),
  ("cdses", "CDSes"),
  ("csv", "CSV"),
  ("dot", "DOT"),
  ("dp", "DP"),
  ("gtr", "GTR"),
  ("json", "JSON"),
  ("mat", "MAT"),
  ("n", "number of"),
  ("nuc", "nucleotide"),
  ("nwk", "Newick"),
  ("opt", "optimization"),
  ("pb", "protobuf"),
  ("pc", "pseudocount"),
  ("pi", "pi"),
  ("tsv", "TSV"),
];

pub fn setting_label(key: &str) -> String {
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
