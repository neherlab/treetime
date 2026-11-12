use crate::runs::record::RunRecord;
use app_commands::command::AppCommand;
use app_output::output_plan::OutputSelection;
use itertools::Itertools;
use std::path::{Path, PathBuf};

const STEM_MAX_CHARS: usize = 60;

const RUN_ID_CHARS: usize = 8;

pub fn download_name(record: &RunRecord, relative: &str) -> String {
  format!("{}.{}", download_stem(record), relative.replace('/', "-"))
}

pub fn package_name(record: &RunRecord) -> String {
  format!("{}.zip", download_stem(record))
}

pub struct FileDownload {
  pub path: PathBuf,
  pub name: String,
}

pub fn output_label(command: AppCommand, kind: Option<OutputSelection>, relative: &str) -> String {
  let Some(kind) = kind else {
    return Path::new(relative)
      .file_name()
      .map_or_else(|| relative.to_owned(), |name| name.to_string_lossy().into_owned());
  };
  match (command, kind) {
    (AppCommand::Timetree, OutputSelection::Nwk) => "Dated tree (Newick, branch lengths in years)".to_owned(),
    (AppCommand::Timetree, OutputSelection::Nexus) => "Dated tree with node annotations (Nexus)".to_owned(),
    (AppCommand::Clock, OutputSelection::Nwk) => "Rerooted tree (Newick)".to_owned(),
    (AppCommand::Optimize, OutputSelection::Nwk) => "Optimized tree (Newick)".to_owned(),
    (AppCommand::Prune, OutputSelection::Nwk) => "Pruned tree (Newick)".to_owned(),
    (_, OutputSelection::RttChart) => format!("Root-to-tip chart ({})", chart_format(relative)),
    (_, kind) => generic_label(kind).to_owned(),
  }
}

pub fn main_outputs(command: AppCommand) -> &'static [OutputSelection] {
  match command {
    AppCommand::Timetree => &[OutputSelection::Auspice, OutputSelection::Nwk, OutputSelection::Nexus],
    AppCommand::Clock => &[OutputSelection::Rtt, OutputSelection::Nwk],
    AppCommand::Ancestral => &[
      OutputSelection::Sequences,
      OutputSelection::Nexus,
      OutputSelection::Auspice,
    ],
    AppCommand::Mugration => &[
      OutputSelection::Traits,
      OutputSelection::Nexus,
      OutputSelection::Auspice,
    ],
    AppCommand::Homoplasy => &[OutputSelection::HomoplasyReport, OutputSelection::HomoplasyStats],
    AppCommand::Optimize | AppCommand::Prune => &[OutputSelection::Nwk],
  }
}

#[rustfmt::skip]
const fn generic_label(kind: OutputSelection) -> &'static str {
  match kind {
    OutputSelection::All                => "All outputs",
    OutputSelection::Nwk                => "Tree (Newick)",
    OutputSelection::Nexus              => "Annotated tree (Nexus)",
    OutputSelection::Auspice            => "Auspice JSON",
    OutputSelection::MatPb              => "Mutation-annotated tree (UShER protobuf)",
    OutputSelection::MatJson            => "Mutation-annotated tree (JSON)",
    OutputSelection::GraphJson          => "Tree graph (JSON)",
    OutputSelection::Dot                => "Tree graph (Graphviz DOT)",
    OutputSelection::NodeData           => "Node data (augur JSON)",
    OutputSelection::Model              => "Substitution models (JSON)",
    OutputSelection::ClockModel         => "Clock model (JSON)",
    OutputSelection::Dates              => "Node dates (TSV)",
    OutputSelection::TraitProbabilities => "Trait probabilities (TSV)",
    OutputSelection::Sequences          => "Ancestral sequences (FASTA)",
    OutputSelection::Translations       => "Ancestral protein sequences (FASTA)",
    OutputSelection::Traits             => "Inferred traits (TSV)",
    OutputSelection::Rtt                => "Root-to-tip table (TSV)",
    OutputSelection::Tracelog           => "Convergence trace (TSV)",
    OutputSelection::Coalescent         => "Coalescent time scale (TSV)",
    OutputSelection::CoalescentJson     => "Coalescent model (JSON)",
    OutputSelection::RttChart           => "Root-to-tip chart",
    OutputSelection::HomoplasyStats     => "Homoplasy statistics (JSON)",
    OutputSelection::HomoplasyReport    => "Homoplasy report (text)",
  }
}

fn chart_format(relative: &str) -> String {
  Path::new(relative).extension().map_or_else(
    || "SVG".to_owned(),
    |extension| extension.to_string_lossy().to_uppercase(),
  )
}

fn download_stem(record: &RunRecord) -> String {
  let title = record
    .title
    .split(|c: char| !c.is_alphanumeric())
    .filter(|word| !word.is_empty())
    .join("-")
    .to_lowercase();
  let title = title
    .chars()
    .take(STEM_MAX_CHARS)
    .collect::<String>()
    .trim_end_matches('-')
    .to_owned();
  let id = record.id.as_str().chars().take(RUN_ID_CHARS).collect::<String>();
  if title.is_empty() {
    format!("treetime-{id}")
  } else {
    format!("{title}-{id}")
  }
}
