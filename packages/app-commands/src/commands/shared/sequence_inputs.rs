use crate::commands::shared::alignment::{PairedAlignment, pair_alignment};
use crate::commands::shared::partition::resolve_run_partitions;
use crate::commands::shared::tree_input::read_input_tree;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::attach::{MissingLeafLimit, complete_alignment_for_leaves};
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::resolve::{PartitionEntry, PartitionPlan};
use treetime::partition::split::SequenceInput;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, SequenceSource};
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node_names::NodeNames;
use treetime_io::fasta::fasta_read_file;
use treetime_io::nwk::{NewickDialect, NwkParse};
use treetime_schema::schema_defaults;

pub(crate) fn read_tree_and_sources(
  tree: &Path,
  newick_dialect: NewickDialect,
  sequences: &SequenceInputs,
  model: &GtrModel,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<TreeAndSources, Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");
  let tree = read_input_tree(tree, newick_dialect, log)?;
  let SequenceSources { alignments, partitions } = sequences.read_partitions(model, &tree.graph, &tree.names, log)?;
  Ok(TreeAndSources {
    tree,
    sources: alignments,
    partitions,
  })
}

pub(crate) fn ancestral_input(
  graph: Graph,
  branch_lengths: EdgeMap<Option<f64>>,
  sequences: SequenceInput,
) -> AncestralInput {
  let edges = branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  AncestralInput {
    graph,
    edges,
    sequences,
  }
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct MissingSequenceArgs {
  /// Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.
  ///
  /// Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may
  /// not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, help_heading = "Input data"))]
  pub allow_missing_sequences: bool,
}

pub fn missing_leaf_limit(allow_missing_sequences: bool) -> MissingLeafLimit {
  MissingLeafLimit {
    ignored: allow_missing_sequences,
  }
}

#[derive(Debug, Clone)]
pub(crate) struct SequenceInputs {
  pub(crate) alignments: Vec<PathBuf>,
  pub(crate) alphabet: Option<AlphabetName>,
  pub(crate) gap_fill: GapFill,
  pub(crate) missing_leaves: MissingLeafLimit,
  pub(crate) partition_entries: Vec<PartitionEntry>,
}

impl SequenceInputs {
  pub(crate) fn read(
    &self,
    graph: &Graph,
    names: &NodeNames,
    log: &dyn LogSink,
  ) -> Result<Vec<AlignmentSource>, Report> {
    let alphabet = Alphabet::new(self.alphabet.unwrap_or_default())?;
    let settings = SourceSettings {
      alphabet: &alphabet,
      gap_fill: Some(self.gap_fill),
      missing_leaves: self.missing_leaves,
    };
    read_sequence_sources(&self.alignments, &settings, graph, names, log)
  }

  pub(crate) fn read_partitions(
    &self,
    model: &GtrModel,
    graph: &Graph,
    names: &NodeNames,
    log: &dyn LogSink,
  ) -> Result<SequenceSources, Report> {
    let alignments = self.read(graph, names, log)?;
    let lengths = alignments.iter().map(|alignment| alignment.source.length).collect_vec();
    let partitions = resolve_run_partitions(&self.partition_entries, model, &lengths, self.alphabet, log)?;
    Ok(SequenceSources { alignments, partitions })
  }
}

pub(crate) fn read_sequence_sources(
  alignments: &[PathBuf],
  settings: &SourceSettings<'_>,
  graph: &Graph,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<Vec<AlignmentSource>, Report> {
  alignments
    .iter()
    .map(|path| read_sequence_source(path, settings, graph, names, log))
    .collect()
}

pub(crate) struct TreeAndSources {
  pub(crate) tree: NwkParse,
  pub(crate) sources: Vec<AlignmentSource>,
  pub(crate) partitions: Vec<PartitionPlan>,
}

pub(crate) struct SequenceSources {
  pub(crate) alignments: Vec<AlignmentSource>,
  pub(crate) partitions: Vec<PartitionPlan>,
}

fn read_sequence_source(
  path: &Path,
  settings: &SourceSettings<'_>,
  graph: &Graph,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<AlignmentSource, Report> {
  let alphabet = settings.alphabet;
  let mut records = fasta_read_file(path, alphabet)?;
  if let Some(gap_fill) = settings.gap_fill {
    for record in &mut records {
      apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
    }
  }
  let n_records = records.len();
  let PairedAlignment { mut sequences, descs } = pair_alignment(records, path, graph, names, log);
  let in_file = || format!("In alignment '{}'", path.display());
  let length = sequences.common_length().wrap_err_with(in_file)?;
  complete_alignment_for_leaves(
    graph,
    &mut sequences.nodes,
    length,
    alphabet,
    settings.missing_leaves,
    &format!("alignment '{}'", path.display()),
    log,
  )
  .wrap_err_with(in_file)?;
  let mask = sequences.mask(length, alphabet);
  Ok(AlignmentSource {
    source: SequenceSource {
      path: path.to_path_buf(),
      alphabet: alphabet.clone(),
      length,
      nodes: sequences.nodes,
      mask,
    },
    descs,
    n_records,
  })
}

pub(crate) struct SourceSettings<'a> {
  pub(crate) alphabet: &'a Alphabet,
  pub(crate) gap_fill: Option<GapFill>,
  pub(crate) missing_leaves: MissingLeafLimit,
}

pub(crate) struct AlignmentSource {
  pub(crate) source: SequenceSource,
  pub(crate) descs: NodeMap<Option<String>>,
  pub(crate) n_records: usize,
}
