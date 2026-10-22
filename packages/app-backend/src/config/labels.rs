use itertools::Itertools;

const LABELS: [(&str, &str); 9] = [
  ("clock_std_dev", "Clock rate std. dev."),
  ("coalescent", "Coalescent prior"),
  ("max_iter", "Iterations"),
  ("method_anc", "Ancestral method"),
  ("model", "Substitution model"),
  ("n_skyline", "Skyline grid points"),
  ("output_selection", "Output files"),
  ("relax", "Relaxed clock (slack, coupling)"),
  ("reroot", "Reroot method"),
];

const WORDS: [(&str, &str); 16] = [
  ("aa", "amino-acid"),
  ("anc", "ancestral"),
  ("cds", "CDS"),
  ("csv", "CSV"),
  ("dot", "DOT"),
  ("gtr", "GTR"),
  ("json", "JSON"),
  ("mat", "MAT"),
  ("n", "number of"),
  ("nuc", "nucleotide"),
  ("nwk", "Newick"),
  ("pb", "protobuf"),
  ("pc", "pseudocount"),
  ("pi", "pi"),
  ("rtt", "root-to-tip"),
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
