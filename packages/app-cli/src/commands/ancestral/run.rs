use crate::commands::ancestral::aa_node_data::{
  read_aa_root_sequences, read_gff3_annotations, template_has_cds_placeholder, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::{TreetimeAncestralArgs, ancestral_params};
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json_with_aa;
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::EdgeMutationCommentProvider;
use app_output::ancestral_result::{AncestralNodeOut, AncestralOutputMaps, AncestralResult, AugurOutputMaps, EdgeOut};
use app_output::ancestral_tree_output::write_ancestral_tree_outputs;
use eyre::Report;
use log::{info, warn};
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::aa::{AaNodeData, reconstruct_aa};
use treetime::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
use treetime::ancestral::mask::create_mask;
use treetime::ancestral::multi::{MarginalPartitionParams, PartitionPlan};
use treetime::ancestral::pipeline::{self, AncestralPartition};
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::progress::ProgressSink;
use treetime::seq::alignment::{EdgeSeqInput, ReconstructionInput, get_common_length, node_seq_inputs};
use treetime::seq::gap_fill::apply_gap_fill;
use treetime::seq::mutation::MutationTrack;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaReader, FastaWriter, read_many_fasta, read_many_fasta_path};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::AlignmentRecord;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub fn run_ancestral_reconstruction(
  args: &TreetimeAncestralArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  validate_aa_args(
    &args.translations,
    &args.cdses,
    &args.annotation,
    &args.aa_root_sequence,
  )?;

  let AncestralInput {
    mut input,
    mask,
    alphabet,
    descs,
    confidences,
  } = read_nwk_fasta(args, cancel, progress)?;
  let names = input.names();
  let branch_lengths = input.branch_lengths();

  let topology_order = args.topology_order.resolve_topology_order(&input.graph, &names, None)?;

  let resolved = args.resolve_outputs()?;
  let output_fasta = if resolved
    .non_tree_outputs
    .contains_key(&OutputSelection::ReconstructedNucFasta)
  {
    let path = &resolved.non_tree_outputs[&OutputSelection::ReconstructedNucFasta];
    Some(FastaWriter::new(create_file_or_stdout(path)?))
  } else {
    None
  };

  let params = ancestral_params(args);

  let result = pipeline::run(&params, &input, alphabet, mask, cancel, progress).map_err(|err| err.into_report())?;

  let aa_fasta_template: Option<String> = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedAaFasta)
    .map(|path| path.to_string_lossy().into_owned());

  let aa_result = if let Some(translations) = &args.translations {
    Some(run_aa_reconstructions(
      args,
      translations,
      aa_fasta_template.as_deref(),
      &input.graph,
      &names,
      &branch_lengths,
      cancel,
      progress,
    )?)
  } else {
    // Prerequisite gating: reconstructed AA FASTA needs --translations. An explicit per-file flag
    // is a hard error; the same output reached via selection or `--output-selection=all` is skipped.
    if aa_fasta_template.is_some() {
      if args.output_reconstructed_aa_fasta.is_some() {
        return make_error!("--output-reconstructed-aa-fasta requires --translations");
      }
      warn!("Skipping reconstructed amino-acid FASTA output: --translations not provided");
    }
    None
  };

  let pipeline::AncestralOutputFull { output, partition } = result;
  let pipeline::AncestralOutput {
    gtr,
    model_name,
    mask,
    emitted_nodes,
  } = output;

  // Stream the reconstructed nucleotide FASTA one record at a time, reading each sequence back off the
  // partition in the walk's emission order. `augur_node_sequence` returns the same flag-aware sequence
  // the reconstruction produced (a posterior draw, a tip echo or imputation, or the MAP state), so the
  // FASTA matches the augur node-data JSON and never holds every sequence in memory at once.
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

  // Gather the per-node/per-edge sequence and mutation values off the pipeline-local partition into
  // plain value maps the output writers consume. This is the only place that reads sequences and
  // mutations from the partition; the tree, node-data, and Newick-comment writers read the maps
  // instead. Node and edge keys stay stable through topology ordering, so gathering before it is
  // bit-identical. The collection runs once and only when a selected tree writer reads sequences, so
  // a Graph-JSON-only, DOT-only, or GTR-only request never expands the sparse sequences.
  let tree_maps = collect_ancestral_tree_maps(&resolved.tree_outputs, || {
    gather_ancestral_output_maps(&input.graph, partition.as_ref())
  })?;
  let augur_maps = if resolved.non_tree_outputs.contains_key(&OutputSelection::AugurNodeData) {
    gather_augur_output_maps_opt(&input.graph, partition.as_ref())?
  } else {
    None
  };

  topology_order.apply(&mut input.graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.9, "");

  // Project the per-node name/confidence and per-edge branch length from the input sidecar maps into
  // the keyed value maps the tree writer and the result consume. The writers read sequences and model
  // metadata from the graph data slot; these maps carry only the name, input-branch-support, and
  // branch-length values. Node and edge keys are stable across topology ordering, so the `names` and
  // `branch_lengths` maps gathered earlier stay valid and are reused for the augur and tree writers.
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

  // Split the AA reconstruction into the node-data result and the CDS annotation map: the annotations
  // come from the CLI's GFF parse and are handed to the augur and tree encoders directly, rather than
  // routed back through the core result.
  let aa_node_data = aa_result.as_ref().map(|(node_data, _)| node_data);
  let empty_aa_annotations = BTreeMap::new();
  let aa_annotations = aa_result
    .as_ref()
    .map_or(&empty_aa_annotations, |(_, annotations)| annotations);

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    if let Some(augur_maps) = &augur_maps {
      write_augur_node_data_json_with_aa(
        &input.graph,
        augur_maps,
        &mask,
        &names,
        aa_node_data,
        aa_annotations,
        path,
      )?;
    }
    info!("Wrote augur node data JSON to {}", path.display());
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, model_name);
        write_gtr_json(&gtr_output, path)?;
      },
      None if args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Use --model=infer or --gtr-iterations.");
      },
      None => warn!("Skipping GTR output: no GTR model was fitted (use --model=infer or --gtr-iterations)"),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    write_tree_for_partition(
      &input.graph,
      &nodes,
      &branch_lengths,
      &tree_maps,
      aa_node_data,
      aa_annotations,
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

/// Read and gap-fill the alignment, parse the tree, complete the alignment so every leaf carries a
/// sequence, compute the alignment mask, and merge everything into the reconstruction input.
///
/// The mask is computed over the completed alignment records (matching the reconstruction's view of
/// the alignment) and returned alongside the merged input and the alphabet the pipeline reconstructs
/// over.
/// The parsed reconstruction input plus the CLI-side sidecars the output writers need but the slim
/// core input drops: per-name leaf descriptions (for the reconstructed FASTA) and per-node
/// input-tree branch support (for the tree and node-data writers), both keyed as the writers consume
/// them.
struct AncestralInput {
  input: ReconstructionInput,
  mask: Vec<bool>,
  alphabet: Alphabet,
  descs: BTreeMap<String, Option<String>>,
  confidences: BTreeMap<GraphNodeKey, Option<f64>>,
}

fn read_nwk_fasta(
  args: &TreetimeAncestralArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<AncestralInput, Report> {
  let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let mut aln = if args.alignment.alignment.is_empty() {
    info!("Reading input fasta from standard input");
    let reader = FastaReader::new(open_stdin()?, &alphabet);
    read_many_fasta(reader)?
  } else {
    read_many_fasta_path(&args.alignment.alignment, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  // Descriptions live only on the input leaf FASTA records and the slim core input drops them, so
  // capture a name-keyed description map here (first record wins on a duplicate name) before the
  // records convert. The reconstructed-FASTA writer rebuilds each node's description by matching the
  // node name back to this map.
  let descs = aln.iter().fold(BTreeMap::new(), |mut descs, record| {
    descs
      .entry(record.seq_name.clone())
      .or_insert_with(|| record.desc.clone());
    descs
  });

  cancel.check()?;
  progress.report("Parsing tree", 0.1, "");
  let parse = nwk_read_file(args.tree())?;

  // The input-tree branch support is a parse-time value the slim core input drops. Capture it keyed
  // by node so the tree and node-data writers can project it; keys stay stable through ancestral
  // reconstruction, which never re-roots.
  let confidences = parse.confidences();

  // Tips absent from the alignment become fully-ambiguous sequences, once, before the mask and the
  // merged input are built, so every leaf's node input carries a sequence and attachment finds one by
  // node key.
  let names = parse.names();
  let aln = aln.into_iter().map(AlignmentRecord::from).collect();
  let aln = complete_alignment_for_leaves(&parse.graph, aln, &alphabet, args.ignore_missing_alns, &names)?;
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);

  let graph = parse.graph;
  let nodes = node_seq_inputs(&graph, &names, aln);
  let edges = parse
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  let input = ReconstructionInput { graph, nodes, edges };
  Ok(AncestralInput {
    input,
    mask,
    alphabet,
    descs,
    confidences,
  })
}

fn write_tree_for_partition(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
  aa_annotations: &BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  resolved: &crate::commands::shared::output::ResolvedOutputs,
  partition: Option<&AncestralPartition>,
) -> Result<(), Report> {
  // Sparse and dense reconstructions annotate Newick/Nexus nodes with their inbound mutations; Fitch
  // parsimony and the partition-less case emit no such comments. The comment provider reads the
  // gathered per-edge mutation map; the partition selects only whether to attach it.
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
    aa_node_data,
    aa_annotations,
    &resolved.tree_outputs,
    &providers,
  )
}

/// Collect the reconstructed sequence and mutation maps once, gated on whether a selected tree writer
/// reads them.
///
/// When only topology-only writers are selected (or no tree writer at all), `gather` is not called and
/// empty maps are returned, so the sparse sequences are never expanded. The single collected result is
/// shared across every selected tree writer.
pub(crate) fn collect_ancestral_tree_maps<F>(
  tree_outputs: &BTreeMap<TreeWriteKind, PathBuf>,
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

/// Whether any selected tree output reads the reconstructed sequence and mutation maps.
///
/// Newick, Nexus, Auspice, and UShER MAT writers read the root sequence and per-edge
/// mutations. The internal Graph JSON dump and the Graphviz DOT writer read only topology and branch
/// lengths, so a selection limited to them needs no sequence collection.
pub(crate) fn tree_outputs_need_sequences(tree_outputs: &BTreeMap<TreeWriteKind, PathBuf>) -> bool {
  tree_outputs
    .keys()
    .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot))
}

/// Gather the root sequence and per-edge nucleotide mutations the tree writers read off the ancestral
/// partition.
pub(crate) fn gather_ancestral_output_maps(
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

/// Gather the augur node-data sequences and substitutions for the partition in the graph data slot, or
/// `None` when no partition exists.
fn gather_augur_output_maps_opt(
  graph: &Graph,
  partition: Option<&AncestralPartition>,
) -> Result<Option<AugurOutputMaps>, Report> {
  let Some(partition) = partition else {
    return Ok(None);
  };
  Ok(Some(gather_augur_output_maps(graph, partition)?))
}

/// Gather the augur node-data sequences and substitutions from one partition.
pub(crate) fn gather_augur_output_maps(
  graph: &Graph,
  partition: &AncestralPartition,
) -> Result<AugurOutputMaps, Report> {
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

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  translations: &str,
  aa_fasta_template: Option<&str>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<(AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>), Report> {
  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let aa_model = ancestral_args.aa_model.resolve();
  let recon_alphabet = Alphabet::new(aa_model.alphabet)?;
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();

  let annotations = read_gff3_annotations(ancestral_args.annotation.as_deref(), &ancestral_args.cdses)?;

  let cdses: Vec<String> = if ancestral_args.cdses.is_empty() {
    annotations.keys().cloned().collect()
  } else {
    ancestral_args.cdses.clone()
  };

  if let Some(aa_seq_template) = aa_fasta_template {
    if cdses.len() > 1 && !template_has_cds_placeholder(aa_seq_template) {
      return make_error!(
        "--output-reconstructed-aa-fasta template needs a CDS placeholder when reconstructing multiple CDSes, \
         otherwise each CDS overwrites the same output file"
      );
    }
  }

  let aa_root_sequences = read_aa_root_sequences(ancestral_args.aa_root_sequence.as_deref(), &cdses, &recon_alphabet)?;

  cancel.check()?;
  progress.report("AA ancestral reconstruction", 0.75, "");

  let params = MarginalPartitionParams {
    dense: ancestral_args.dense,
    include_leaves: ancestral_args.include_leaves || ancestral_args.reconstruct_tip_states,
    impute_missing_data: ancestral_args.impute_missing_data || ancestral_args.reconstruct_tip_states,
    sample_from_profile: ancestral_args.sample_from_profile,
    seed: ancestral_args.seed,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  // Read and sanitize each CDS translation FASTA, then build its reconstruction plan in CDS order.
  // Sanitizing folds out-of-alphabet amino-acid characters (e.g. stop '*') into the reconstruction
  // alphabet's unknown state, and gap-fill applies the same overhang policy as the nucleotide path.
  // The core driver reconstructs the plans in this order.
  let mut plans = Vec::with_capacity(cdses.len());
  for cds in &cdses {
    let path = translation_path(translations, cds);
    let mut sequences = read_many_fasta_path(&[&path], &read_alphabet)?;
    let mut sanitized = 0_usize;
    for record in &mut sequences {
      let (seq, changed) = sanitize_to_alphabet(&record.seq, &recon_alphabet);
      record.seq = seq;
      sanitized += changed;
      apply_gap_fill(
        &mut record.seq,
        gap_fill_mode,
        recon_alphabet.gap(),
        recon_alphabet.unknown(),
      );
    }
    if sanitized > 0 {
      warn!(
        "CDS '{cds}': mapped {sanitized} out-of-alphabet amino-acid characters (e.g. stop '*') to '{}'.",
        char::from(recon_alphabet.unknown())
      );
    }

    plans.push(PartitionPlan {
      name: cds.clone(),
      alphabet: recon_alphabet.clone(),
      gtr_model: aa_model.gtr_model,
      sequences: sequences.into_iter().map(AlignmentRecord::from).collect(),
      annotation: annotations.get(cds).cloned(),
      reference_override: aa_root_sequences.get(cds).cloned(),
    });
  }

  // Per-CDS reconstructed amino-acid FASTA sink: the core driver streams every node's sequence while
  // its CDS partition is resident, and the sink opens a fresh output file when the CDS changes (the
  // driver finishes one CDS before starting the next), writing each node in tree order.
  let seq_sink: Option<Box<dyn SeqSink>> = aa_fasta_template
    .map(|template| -> Box<dyn SeqSink> { Box::new(AaFastaSink::new(template.to_owned(), names.clone())) });

  // The CDS annotation map the augur and tree encoders need, kept CLI-side: the reconstructed subset
  // of the GFF-parsed annotations (one entry per reconstructed CDS that has an annotation). The core
  // result no longer carries it.
  let cds_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry> = cdses
    .iter()
    .filter_map(|cds| annotations.get(cds).map(|entry| (cds.clone(), entry.clone())))
    .collect();

  let node_data = reconstruct_aa(graph, names, branch_lengths, &params, plans, seq_sink)?;
  Ok((node_data, cds_annotations))
}

/// Per-CDS reconstructed amino-acid FASTA sink.
///
/// Opens a fresh output file when the emitted CDS changes and writes each node's sequence in tree
/// order. Node names come from the parsed tree, with a `node_{key}` fallback for unnamed internal nodes
/// (matching the augur node-data naming); amino-acid records carry no description. The amino-acid
/// reconstruction does not change the topology, so the names captured at construction stay valid and
/// `on_topology` needs no per-topology resolution.
struct AaFastaSink {
  template: String,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  open: Option<(String, FastaWriter)>,
}

impl AaFastaSink {
  fn new(template: String, names: BTreeMap<GraphNodeKey, Option<String>>) -> Self {
    Self {
      template,
      names,
      open: None,
    }
  }
}

impl SeqSink for AaFastaSink {
  fn on_topology(&mut self, _graph: &Graph) -> Result<(), Report> {
    Ok(())
  }

  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    let SeqTrack::Aa(cds) = item.track else {
      return treetime_utils::make_internal_error!("Amino-acid reconstructed FASTA sink received a nucleotide track");
    };
    if self.open.as_ref().map(|(name, _)| name.as_str()) != Some(cds) {
      let path = translation_path(&self.template, cds);
      self.open = Some((cds.to_owned(), FastaWriter::new(create_file_or_stdout(path)?)));
    }
    let name = self.names[&item.key]
      .as_deref()
      .map_or_else(|| format!("node_{}", item.key.0), str::to_owned);
    self
      .open
      .as_mut()
      .expect("writer opened above")
      .1
      .write(&name, &None, item.seq)
  }
}
