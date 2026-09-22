use crate::commands::support::{default_output_plan, default_topology_order};
use app_output::EdgeMutationCommentProvider;
use app_output::ancestral_result::{AncestralNodeOut, AncestralOutputMaps, AncestralResult, AugurOutputMaps, EdgeOut};
use app_output::ancestral_tree_output::write_ancestral_tree_outputs;
use app_output::augur_node_data_ancestral::write_augur_node_data_json_with_aa;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use eyre::Report;
use log::{info, warn};
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::attach::complete_alignment_for_leaves;
use treetime::ancestral::mask::create_mask;
use treetime::ancestral::params::MethodAncestral;
use treetime::ancestral::pipeline::{self, AncestralParams, AncestralPartition};
use treetime::ancestral::sample::SampleMode;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::progress::ProgressSink;
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime::seq::mutation::MutationTrack;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaReader, FastaWriter, read_many_fasta, read_many_fasta_path};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::{CommentProviders, nwk_read_file};
use treetime_primitives::AlignmentRecord;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use utoipa::ToSchema;

/// Ancestral reconstruction request (openapi subset).
#[derive(Debug, SmartDefault, Deserialize, ToSchema)]
#[serde(default)]
pub struct AncestralArgs {
  input_fastas: Vec<String>,
  aln: Option<String>,
  vcf_reference: Option<String>,
  tree: String,
  #[schema(value_type = Option<String>)]
  alphabet: Option<AlphabetName>,
  #[default(GtrModelName::Infer)]
  #[schema(value_type = String)]
  model_name: GtrModelName,
  gtr_params: Vec<String>,
  #[default(MethodAncestral::default())]
  #[schema(value_type = String)]
  method_anc: MethodAncestral,
  dense: Option<bool>,
  aa: bool,
  #[default(GapFill::default())]
  #[schema(value_type = String)]
  gap_fill: GapFill,
  keep_overhangs: bool,
  zero_based: bool,
  include_leaves: bool,
  impute_missing_data: bool,
  reconstruct_tip_states: bool,
  report_ambiguous: bool,
  outdir: String,
  gtr_iterations: usize,
  site_specific_gtr: bool,
  seed: Option<u64>,
}

impl AncestralArgs {
  fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }

  fn params(&self) -> AncestralParams {
    AncestralParams {
      method: self.method_anc,
      model: self.model_name,
      dense: self.dense,
      include_leaves: self.include_leaves || self.reconstruct_tip_states,
      impute_missing_data: self.impute_missing_data || self.reconstruct_tip_states,
      gtr_iterations: self.gtr_iterations,
      site_specific_gtr: self.site_specific_gtr,
      seed: self.seed,
      sample_from_profile: SampleMode::default(),
      ignore_missing_alns: false,
    }
  }
}

pub(crate) fn run_ancestral(
  args: &AncestralArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  let AncestralReadInputs {
    mut input,
    mask,
    descs,
    confidences,
  } = read_nwk_fasta(args, cancel, progress)?;
  let names = input.names();
  let branch_lengths = input.branch_lengths();

  let resolved = default_output_plan(CommandKind::Ancestral, std::path::Path::new(&args.outdir))?;
  let output_fasta = if resolved
    .non_tree_outputs
    .contains_key(&OutputSelection::ReconstructedNucFasta)
  {
    let path = &resolved.non_tree_outputs[&OutputSelection::ReconstructedNucFasta];
    Some(FastaWriter::new(create_file_or_stdout(path)?))
  } else {
    None
  };

  let params = args.params();
  let alphabet = Alphabet::new(args.alphabet.unwrap_or_default())?;

  let result = pipeline::run(&params, &input, alphabet, mask, cancel, progress).map_err(|err| err.into_report())?;

  let pipeline::AncestralOutputFull { output, partition } = result;
  let pipeline::AncestralOutput {
    gtr,
    model_name,
    mask,
    emitted_nodes,
  } = output;

  if let Some(mut writer) = output_fasta {
    if let Some(partition) = partition.as_ref() {
      for &key in &emitted_nodes {
        let node = &input.nodes[&key];
        let desc = node.name.as_deref().and_then(|name| descs.get(name)).cloned().flatten();
        writer.write(
          node.name.as_deref().unwrap_or(""),
          &desc,
          &partition.augur_node_sequence(key),
        )?;
      }
    }
  }

  let tree_maps = collect_ancestral_tree_maps(&resolved.tree_outputs, || {
    gather_ancestral_output_maps(&input.graph, partition.as_ref())
  })?;
  let augur_maps = if resolved.non_tree_outputs.contains_key(&OutputSelection::AugurNodeData) {
    gather_augur_output_maps_opt(&input.graph, partition.as_ref())?
  } else {
    None
  };

  default_topology_order().apply(&mut input.graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.9, "");

  let nodes: BTreeMap<GraphNodeKey, AncestralNodeOut> = input
    .nodes
    .iter()
    .map(|(key, node)| {
      (
        *key,
        AncestralNodeOut {
          name: node.name.clone(),
          confidence: confidences.get(key).copied().flatten(),
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = input
    .edges
    .iter()
    .map(|(key, edge)| {
      (
        *key,
        EdgeOut {
          branch_length: edge.branch_length,
        },
      )
    })
    .collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    if let Some(augur_maps) = &augur_maps {
      write_augur_node_data_json_with_aa(&input.graph, augur_maps, &mask, &names, None, &BTreeMap::new(), path)?;
    }
    info!("Wrote augur node data JSON to {}", path.display());
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, model_name);
        write_gtr_json(&gtr_output, path)?;
      },
      None => warn!("Skipping GTR output: no GTR model was fitted (use model=infer or gtr_iterations)"),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    write_tree_for_partition(
      &input.graph,
      &nodes,
      &branch_lengths,
      &tree_maps,
      &resolved,
      partition.as_ref(),
    )?;
  }

  progress.report("Done", 1.0, "");
  Ok(AncestralResult {
    graph: input.graph,
    nodes,
    edges,
  })
}

struct AncestralReadInputs {
  input: AncestralInput,
  mask: Vec<bool>,
  descs: BTreeMap<String, Option<String>>,
  confidences: BTreeMap<GraphNodeKey, Option<f64>>,
}

fn read_nwk_fasta(
  args: &AncestralArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<AncestralReadInputs, Report> {
  let gap_fill_mode = args.effective_gap_fill();
  let alphabet = Alphabet::new(args.alphabet.unwrap_or_default())?;

  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let mut aln = if args.input_fastas.is_empty() {
    info!("Reading input fasta from standard input");
    let reader = FastaReader::new(open_stdin()?, &alphabet);
    read_many_fasta(reader)?
  } else {
    let paths: Vec<std::path::PathBuf> = args.input_fastas.iter().map(std::path::PathBuf::from).collect();
    read_many_fasta_path(&paths, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  let descs = aln.iter().fold(BTreeMap::new(), |mut descs, record| {
    descs
      .entry(record.seq_name.clone())
      .or_insert_with(|| record.desc.clone());
    descs
  });

  cancel.check()?;
  progress.report("Parsing tree", 0.1, "");
  let parse = nwk_read_file(std::path::Path::new(&args.tree))?;
  let confidences = parse.confidences();

  let names = parse.names();
  let aln = aln.into_iter().map(AlignmentRecord::from).collect();
  let aln = complete_alignment_for_leaves(&parse.graph, aln, &alphabet, false, &names)?;
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);

  let graph = parse.graph;
  let nodes = node_seq_inputs(&graph, &names, aln);
  let edges = parse
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  let input = AncestralInput { graph, nodes, edges };
  Ok(AncestralReadInputs {
    input,
    mask,
    descs,
    confidences,
  })
}

fn write_tree_for_partition(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  resolved: &ResolvedOutputs,
  partition: Option<&AncestralPartition>,
) -> Result<(), Report> {
  let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, graph);
  let providers = match partition {
    Some(AncestralPartition::Sparse(_) | AncestralPartition::Dense(_)) => CommentProviders::new().with(&provider),
    Some(AncestralPartition::Fitch(_)) | None => CommentProviders::new(),
  };
  write_ancestral_tree_outputs(
    graph,
    nodes,
    branch_lengths,
    maps,
    None,
    &BTreeMap::new(),
    &resolved.tree_outputs,
    &providers,
  )
}

fn collect_ancestral_tree_maps<F>(
  tree_outputs: &BTreeMap<TreeWriteKind, std::path::PathBuf>,
  gather: F,
) -> Result<AncestralOutputMaps, Report>
where
  F: FnOnce() -> Result<AncestralOutputMaps, Report>,
{
  if tree_outputs_need_sequences(tree_outputs) {
    gather()
  } else {
    Ok(AncestralOutputMaps::default())
  }
}

fn tree_outputs_need_sequences(tree_outputs: &BTreeMap<TreeWriteKind, std::path::PathBuf>) -> bool {
  tree_outputs
    .keys()
    .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot))
}

fn gather_ancestral_output_maps(
  graph: &Graph,
  partition: Option<&AncestralPartition>,
) -> Result<AncestralOutputMaps, Report> {
  let Some(partition) = partition else {
    return Ok(AncestralOutputMaps::default());
  };
  let root_sequence = Some(partition.root_sequence(graph)?);
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((key, partition.edge_mutations(graph, key, &MutationTrack::Nucleotide)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(AncestralOutputMaps {
    root_sequence,
    edge_mutations,
  })
}

fn gather_augur_output_maps_opt(
  graph: &Graph,
  partition: Option<&AncestralPartition>,
) -> Result<Option<AugurOutputMaps>, Report> {
  let Some(partition) = partition else {
    return Ok(None);
  };
  Ok(Some(gather_augur_output_maps(graph, partition)?))
}

fn gather_augur_output_maps(graph: &Graph, partition: &AncestralPartition) -> Result<AugurOutputMaps, Report> {
  let sequence_length = partition.sequence_length();
  let ambiguous_char = partition.ambiguous_char();
  let root_sequence = partition.augur_root_sequence(graph)?;
  let node_sequences = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, partition.augur_node_sequence(key))
    })
    .collect();
  let edge_subs = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((key, partition.edge_subs(graph, key)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(AugurOutputMaps {
    root_sequence,
    node_sequences,
    edge_subs,
    sequence_length,
    ambiguous_char,
  })
}
