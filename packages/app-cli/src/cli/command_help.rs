use itertools::Itertools;

pub(crate) const HELP_HINT_OUTPUTS: &str =
  "Use --help to see all options, including the per-file output flags (--output-*)";

pub(crate) const HELP_HINT: &str = "Use --help to see all options";

pub(crate) const TIMETREE_USAGE: &str = "treetime timetree --tree <FILE> --metadata <FILE> [(--alignment <FILE>)... | --sequence-length <N>] (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime timetree --config <FILE> [OPTIONS]";

pub(crate) const TIMETREE_EXAMPLES: &[&str] = &[
  "treetime timetree -t tree.nwk -m metadata.tsv -a aln.fasta.xz -O results",
  "treetime timetree -t tree.nwk -m metadata.tsv --sequence-length 10700 --coalescent opt --confidence -O results",
  "treetime timetree --config timetree.yaml -O results",
];

pub(crate) const ANCESTRAL_USAGE: &str =
  "treetime ancestral --tree <FILE> (--alignment <FILE>)... (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime ancestral --config <FILE> [OPTIONS]";

pub(crate) const ANCESTRAL_EXAMPLES: &[&str] = &[
  "treetime ancestral -t tree.nwk -a aln.fasta.xz -O results",
  "treetime ancestral -t tree.nwk -a ha.fasta.xz -a na.fasta.xz --reconstruct-tip-states --output-sequences 'sequences-{alignment}.fasta'",
  "treetime ancestral -t tree.nwk -a aln.fasta.xz --annotation genome.gff3 --translations 'translations/{cds}.fasta' -O results",
];

pub(crate) const CLOCK_USAGE: &str =
  "treetime clock --tree <FILE> --metadata <FILE> (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime clock --config <FILE> [OPTIONS]";

pub(crate) const CLOCK_EXAMPLES: &[&str] = &[
  "treetime clock -t tree.nwk -m metadata.tsv -O results",
  "treetime clock -t tree.nwk -m metadata.tsv --keep-root --output-rtt rtt.tsv",
];

pub(crate) const MUGRATION_USAGE: &str =
  "treetime mugration --tree <FILE> --metadata <FILE> --attribute <COLUMN> (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime mugration --config <FILE> [OPTIONS]";

pub(crate) const MUGRATION_EXAMPLES: &[&str] = &[
  "treetime mugration -t tree.nwk -m metadata.tsv --attribute country -O results",
  "treetime mugration -t tree.nwk -m metadata.tsv --attribute country --weights weights.tsv --confidence -O results",
];

pub(crate) const HOMOPLASY_USAGE: &str =
  "treetime homoplasy --tree <FILE> --alignment <FILE> (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime homoplasy --config <FILE> [OPTIONS]";

pub(crate) const HOMOPLASY_EXAMPLES: &[&str] = &[
  "treetime homoplasy -t tree.nwk -a aln.fasta.xz -O results",
  "treetime homoplasy -t tree.nwk -a aln.fasta.xz --drms drms.tsv --detailed --output-homoplasy-report report.txt",
];

pub(crate) const OPTIMIZE_USAGE: &str =
  "treetime optimize --tree <FILE> (--alignment <FILE>)... (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime optimize --config <FILE> [OPTIONS]";

pub(crate) const OPTIMIZE_EXAMPLES: &[&str] = &[
  "treetime optimize -t tree.nwk -a aln.fasta.xz -O results",
  "treetime optimize -t tree.nwk -a aln.fasta.xz --model hky85 --reroot min-dev --output-tree-nwk optimized.nwk",
];

pub(crate) const PRUNE_USAGE: &str =
  "treetime prune --tree <FILE> [--alignment <FILE>]... (-O <DIR> | --output-* <FILE>) [OPTIONS]
       treetime prune --config <FILE> [OPTIONS]";

pub(crate) const PRUNE_EXAMPLES: &[&str] = &[
  "treetime prune -t tree.nwk --collapse-short 1e-6 --output-tree-nwk pruned.nwk",
  "treetime prune -t tree.nwk -a aln.fasta.xz --collapse-empty --merge-shared-mutations -O results",
  "treetime prune -t tree.nwk --remove-file outliers.txt --output-tree-nwk pruned.nwk",
];

pub(crate) const PIPELINE_USAGE: &str = "treetime pipeline --config <FILE> [-O <DIR>] [OPTIONS]";

pub(crate) const PIPELINE_EXAMPLES: &[&str] = &[
  "treetime pipeline --config pipeline.yaml -O results",
  "treetime pipeline --config pipeline.yaml --steps clock,timetree -O results",
  "treetime pipeline --config pipeline.yaml --dry-run",
];

pub(crate) const SCHEMA_WRITE_USAGE: &str = "treetime schema write [--for <TARGET>] [-o <PATH>] [OPTIONS]";

pub(crate) const SCHEMA_WRITE_EXAMPLES: &[&str] = &[
  "treetime schema write --for timetree -o timetree.schema.json",
  "treetime schema write --for all -o schemas",
];

pub(crate) const COMPLETIONS_USAGE: &str = "treetime completions [SHELL] [OPTIONS]";

pub(crate) const COMPLETIONS_EXAMPLES: &[&str] = &[
  "treetime completions bash > ~/.local/share/bash-completion/completions/treetime",
  "treetime completions zsh > ~/.zfunc/_treetime",
];

pub(crate) const HELP_MARKDOWN_USAGE: &str = "treetime help-markdown [OPTIONS] > <FILE>";

pub(crate) const HELP_MARKDOWN_EXAMPLES: &[&str] =
  &["treetime help-markdown > reference.md", "treetime help-markdown | less"];

pub(crate) fn after_help(examples: &[&str], hint: &str) -> String {
  format!("{}\n\n{hint}", after_long_help(examples))
}

pub(crate) fn after_long_help(examples: &[&str]) -> String {
  let examples = examples.iter().map(|example| format!("  {example}")).join("\n");
  format!("Examples:\n{examples}")
}
