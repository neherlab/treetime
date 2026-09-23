use crate::commands::support::{default_output_plan, default_topology_order};
use app_output::EdgeMutationCommentProvider;
use app_output::augur_node_data_optimize::write_augur_node_data_json;
use app_output::optimize_result::{EdgeOut, OptimizeNodeOut, OptimizeOutputMaps, OptimizeResult};
use app_output::optimize_tree_output::write_optimize_tree_outputs;
use app_output::output_plan::{CommandKind, OutputSelection};
use eyre::Report;
use log::info;
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use treetime::cancel::Cancel;
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::progress::ProgressSink;
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::nwk::{CommentProviders, nwk_read_file};
use treetime_primitives::{AlignmentRecord, Seq};

pub fn run_optimize(
  args: &OptimizeArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<OptimizeResult, Report> {
  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let alphabet = Alphabet::new(args.alphabet.unwrap_or_default())?;
  let gap_fill = args.effective_gap_fill();
  let paths: Vec<std::path::PathBuf> = args.input_fastas.iter().map(std::path::PathBuf::from).collect();
  let mut aln = read_many_fasta_path(&paths, &alphabet)?;
  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
  }
  let nwk_parsed = nwk_read_file(Path::new(&args.tree))?;
  let confidences = nwk_parsed.confidences();
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;

  let resolved = default_output_plan(CommandKind::Optimize, Path::new(&args.outdir))?;

  let params = OptimizeParams {
    model: args.model_name,
    dense: args.dense,
    max_iter: args.max_iter,
    dp: args.dp,
    damping: args.damping,
    opt_method: args.opt_method,
    initial_guess: args.branch_length_initial_guess,
    no_indels: args.no_indels,
    reroot_spec: args.reroot_spec(),
    topology_ops: TopologyOps::default(),
  };

  let input = OptimizeInput {
    graph,
    alphabet,
    sequences: aln.into_iter().map(AlignmentRecord::from).collect(),
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, cancel, progress).map_err(|err| err.into_report())?;
  let pipeline::OptimizeOutput {
    mut graph,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
    branch_lengths,
    names,
  } = output;

  let maps = gather_optimize_output_maps(&graph, &sparse_partitions, &dense_partitions)?;

  default_topology_order().apply(&mut graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.9, "");

  let nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut> = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (
        key,
        OptimizeNodeOut {
          name: names[&key].clone(),
          confidence: confidences.get(&key).copied().flatten(),
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = branch_lengths
    .iter()
    .map(|(&key, &branch_length)| (key, EdgeOut { branch_length }))
    .collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::builder().gtr(&gtr).model_name(model_name).build();
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    if maps.root_sequence.is_some() {
      let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, &graph);
      let providers = CommentProviders::new().with(&provider);
      write_optimize_tree_outputs(
        &graph,
        &nodes,
        &branch_lengths,
        &maps,
        &resolved.tree_outputs,
        &providers,
      )?;
    } else {
      write_optimize_tree_outputs(
        &graph,
        &nodes,
        &branch_lengths,
        &maps,
        &resolved.tree_outputs,
        &CommentProviders::new(),
      )?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let alignment = args.input_fastas.first().map(Path::new);
    write_augur_node_data_json(
      &graph,
      &nodes,
      &branch_lengths,
      alignment,
      Some(Path::new(&args.tree)),
      None,
      path,
    )?;
    info!("Wrote augur node data JSON to {path}", path = path.display());
  }

  progress.report("Done", 1.0, "");

  Ok(OptimizeResult { graph, nodes, edges })
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct OptimizeArgs {
  pub input_fastas: Vec<String>,
  pub tree: String,
  pub alphabet: Option<AlphabetName>,
  #[default(GtrModelName::Infer)]
  pub model_name: GtrModelName,
  pub dense: Option<bool>,
  pub outdir: String,
  #[default = 10]
  pub max_iter: usize,
  #[default = 0.1]
  pub dp: f64,
  #[default = 0.75]
  pub damping: f64,
  #[default(InitialGuessMode::Auto)]
  pub branch_length_initial_guess: InitialGuessMode,
  #[default(BranchOptMethod::default())]
  pub opt_method: BranchOptMethod,
  pub no_indels: bool,
  pub reroot: Option<OptimizeRerootMethod>,
  pub reroot_tips: Vec<String>,
  pub keep_root: bool,
  #[default(GapFill::default())]
  pub gap_fill: GapFill,
  pub keep_overhangs: bool,
}

impl OptimizeArgs {
  fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }

  fn reroot_spec(&self) -> Option<RerootSpec> {
    if self.keep_root {
      return None;
    }
    if let Some(method) = self.reroot {
      return Some(RerootSpec::Method(RerootMethod::from(method)));
    }
    if !self.reroot_tips.is_empty() {
      return Some(RerootSpec::Tips(self.reroot_tips.clone()));
    }
    None
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum OptimizeRerootMethod {
  MinDev,
}

impl From<OptimizeRerootMethod> for RerootMethod {
  fn from(m: OptimizeRerootMethod) -> Self {
    match m {
      OptimizeRerootMethod::MinDev => RerootMethod::MinDev,
    }
  }
}

fn gather_optimize_output_maps(
  graph: &Graph,
  sparse_partitions: &[SparseReconstruction],
  dense_partitions: &[DenseReconstruction],
) -> Result<OptimizeOutputMaps, Report> {
  if let Some(family) = dense_partitions.first() {
    gather_optimize_partition_maps(
      graph,
      family.root_sequence(graph)?,
      |key| family.edge_mutations(graph, key, &MutationTrack::Nucleotide),
      |key| family.edge_subs(graph, key),
    )
  } else if let Some(family) = sparse_partitions.first() {
    gather_optimize_partition_maps(
      graph,
      family.root_sequence(graph)?,
      |key| family.edge_mutations(key, &MutationTrack::Nucleotide),
      |key| family.edge_subs(key),
    )
  } else {
    Ok(OptimizeOutputMaps::default())
  }
}

fn gather_optimize_partition_maps(
  graph: &Graph,
  root_sequence: Seq,
  edge_mutations: impl Fn(GraphEdgeKey) -> Result<Vec<Mutation>, Report>,
  edge_subs: impl Fn(GraphEdgeKey) -> Result<Vec<Sub>, Report>,
) -> Result<OptimizeOutputMaps, Report> {
  let root_sequence = Some(root_sequence);
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((key, edge_mutations(key)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let edge_subs = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((key, edge_subs(key)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(OptimizeOutputMaps {
    root_sequence,
    edge_mutations,
    edge_subs,
  })
}
