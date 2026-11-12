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

const OPTION_LABELS: [(&str, &str); 26] = [
  ("aa", "Amino acid"),
  ("aa-no-stop", "Amino acid without stop"),
  ("argmax", "Most likely state"),
  ("auto", "Automatic"),
  ("beast", "BEAST"),
  ("classic", "Classic Newick"),
  ("enewick", "Extended Newick"),
  ("f81", "F81"),
  ("hky85", "HKY85"),
  ("infer", "Infer GTR"),
  ("jc69", "JC69"),
  ("jtt92", "JTT92"),
  ("k80", "K80"),
  ("keep", "Input order"),
  ("min-dev", "Minimum deviation"),
  ("mrbayes", "MrBayes"),
  ("nhx", "NHX"),
  ("nuc", "Nucleotide"),
  ("oldest", "Oldest tip"),
  ("only-terminal", "Only terminal gaps"),
  ("opt", "Optimized constant"),
  ("rich", "Rich Newick"),
  ("root", "Sample at the root"),
  ("t92", "T92"),
  ("tn93", "TN93"),
  ("tracelog", "Trace log"),
];

const OPTION_PART_SEPARATOR: char = ',';

const REVERSED_SUFFIX: &str = "-reverse";

pub fn setting_label(key: &str) -> String {
  LABELS
    .iter()
    .find(|(known, _)| *known == key)
    .map_or_else(|| sentence_of_words(key.split('_')), |(_, label)| (*label).to_owned())
}

pub fn option_label(value: &str) -> String {
  value
    .split(OPTION_PART_SEPARATOR)
    .enumerate()
    .map(|(index, part)| {
      let label = option_part_label(part);
      if index == 0 { label } else { in_sentence(&label) }
    })
    .join(", ")
}

fn option_part_label(part: &str) -> String {
  if let Some(base) = part.strip_suffix(REVERSED_SUFFIX) {
    return format!("{}, reversed", option_part_label(base));
  }
  OPTION_LABELS.iter().find(|(known, _)| *known == part).map_or_else(
    || sentence_of_words(part.split(['-', '_'])),
    |(_, label)| (*label).to_owned(),
  )
}

fn sentence_of_words<'a>(words: impl Iterator<Item = &'a str>) -> String {
  let words = words
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

fn in_sentence(label: &str) -> String {
  let mut chars = label.chars();
  match chars.next() {
    Some(first) if !chars.as_str().contains(char::is_uppercase) => first.to_lowercase().chain(chars).collect(),
    _ => label.to_owned(),
  }
}
