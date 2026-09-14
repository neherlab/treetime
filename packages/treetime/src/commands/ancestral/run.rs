use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
use crate::ancestral::mask::create_mask;
use crate::ancestral::multi::{MarginalPartitionParams, PartitionPlan, reconstruct_marginal_partition};
use crate::ancestral::pipeline::{self, AncestralParams, AncestralPartition};
use crate::commands::ancestral::aa_node_data::{
  AaNodeData, annotation_cds_nuc_length, collect_aa_cds_node_data, read_aa_root_sequences, read_gff3_annotations,
  template_has_cds_placeholder, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json_with_aa;
use crate::commands::ancestral::result::{
  AncestralNodeOut, AncestralOutputMaps, AncestralResult, AugurOutputMaps, EdgeOut,
};
use crate::commands::ancestral::tree_output::write_ancestral_tree_outputs;
use crate::commands::shared::mutation_comment::EdgeMutationCommentProvider;
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::progress::ProgressSink;
use crate::seq::alignment::get_common_length;
use crate::seq::gap_fill::apply_gap_fill;
use crate::seq::mutation::MutationTrack;
use eyre::Report;
use log::{info, warn};
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaReader, FastaWriter, read_many_fasta, read_many_fasta_path};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::{NwkFastaInput, nwk_read_file};
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use treetime_utils::sync::random::get_random_number_generator;

pub fn run_ancestral_reconstruction(
  args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  validate_aa_args(
    &args.translations,
    &args.cdses,
    &args.annotation,
    &args.aa_root_sequence,
  )?;

  let (mut input, mask, alphabet) = read_nwk_fasta(args, progress)?;
  let names = input.names();
  let branch_lengths = input.branch_lengths();

  let topology_order = args.topology_order.resolve_topology_order(&input.graph, &names, None)?;

  let resolved = args.resolve_outputs()?;
  let mut output_fasta = if resolved
    .non_tree_outputs
    .contains_key(&OutputSelection::ReconstructedNucFasta)
  {
    let path = &resolved.non_tree_outputs[&OutputSelection::ReconstructedNucFasta];
    Some(FastaWriter::new(create_file_or_stdout(path)?))
  } else {
    None
  };

  let params = AncestralParams::new(args);

  let result = pipeline::run(
    &params,
    &input,
    alphabet,
    mask,
    |key, seq| {
      if let Some(ref mut writer) = output_fasta {
        let node = &input.nodes[&key];
        writer.write(node.name.as_deref().unwrap_or(""), &node.desc, seq)
      } else {
        Ok(())
      }
    },
    progress,
  )?;

  let aa_fasta_template: Option<String> = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedAaFasta)
    .map(|path| path.to_string_lossy().into_owned());

  let aa_node_data = if let Some(translations) = &args.translations {
    Some(run_aa_reconstructions(
      args,
      translations,
      aa_fasta_template.as_deref(),
      &input.graph,
      &names,
      &branch_lengths,
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
    gtr, model_name, mask, ..
  } = output;

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

  // Gather the per-node name/confidence and per-edge branch length off the ordered tree into keyed
  // value maps the output writers consume. The writers still read sequences and model metadata from
  // the graph data slot; these maps carry the name, input-branch-support, and branch-length values.
  let nodes: BTreeMap<GraphNodeKey, AncestralNodeOut> = input
    .graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      let node_input = &input.nodes[&key];
      (
        key,
        AncestralNodeOut {
          name: node_input.name.clone(),
          confidence: node_input.confidence,
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = input
    .graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      (
        key,
        EdgeOut {
          branch_length: input.edges[&key].branch_length,
        },
      )
    })
    .collect();
  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.branch_length)).collect();
  let node_names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    if let Some(augur_maps) = &augur_maps {
      write_augur_node_data_json_with_aa(&input.graph, augur_maps, &mask, &node_names, aa_node_data.as_ref(), path)?;
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
      aa_node_data.as_ref(),
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
fn read_nwk_fasta(
  args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<(NwkFastaInput, Vec<bool>, Alphabet), Report> {
  let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;

  progress.check_cancelled()?;
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

  progress.check_cancelled()?;
  progress.report("Parsing tree", 0.1, "");
  let parse = nwk_read_file(args.tree())?;

  // Tips absent from the alignment become fully-ambiguous sequences, once, before the mask and the
  // merged input are built, so every leaf's node input carries a sequence and attachment finds one by
  // node key.
  let names = parse.names();
  let aln = complete_alignment_for_leaves(&parse.graph, aln, &alphabet, args.ignore_missing_alns, &names)?;
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);

  let input = NwkFastaInput::from_parse_and_aln(parse, aln);
  Ok((input, mask, alphabet))
}

fn write_tree_for_partition(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
  resolved: &crate::commands::shared::output::ResolvedOutputs,
  partition: Option<&AncestralPartition>,
) -> Result<(), Report> {
  // Sparse and dense reconstructions annotate Newick/Nexus nodes with their inbound mutations; Fitch
  // parsimony and the partition-less case emit no such comments. The comment provider reads the
  // gathered per-edge mutation map; the partition selects only whether to attach it.
  match partition {
    Some(AncestralPartition::Sparse(_) | AncestralPartition::Dense(_)) => {
      let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, graph);
      let providers = CommentProviders::new().with(&provider);
      write_ancestral_tree_outputs(
        graph,
        nodes,
        branch_lengths,
        maps,
        aa_node_data,
        &resolved.tree_outputs,
        &providers,
      )?;
    },
    Some(AncestralPartition::Fitch(_)) | None => {
      write_ancestral_tree_outputs(
        graph,
        nodes,
        branch_lengths,
        maps,
        aa_node_data,
        &resolved.tree_outputs,
        &CommentProviders::new(),
      )?;
    },
  }
  Ok(())
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
/// Newick, Nexus, Auspice, PhyloXML, and UShER MAT writers read the per-node sequences and per-edge
/// mutations. The internal Graph JSON dump and the Graphviz DOT writer read only topology and branch
/// lengths, so a selection limited to them needs no sequence collection.
pub(crate) fn tree_outputs_need_sequences(tree_outputs: &BTreeMap<TreeWriteKind, PathBuf>) -> bool {
  tree_outputs
    .keys()
    .any(|kind| !matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot))
}

/// Gather the per-node nucleotide sequences, root sequence, and per-edge nucleotide mutations the tree
/// writers read off the ancestral partition.
pub(crate) fn gather_ancestral_output_maps(
  graph: &Graph,
  partition: Option<&AncestralPartition>,
) -> Result<AncestralOutputMaps, Report> {
  let Some(partition) = partition else {
    return Ok(AncestralOutputMaps::default());
  };
  let root_sequence = Some(partition.root_sequence(graph)?);
  let node_sequences = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.node_sequence(key))
    })
    .collect();
  let edge_mutations = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      Ok((key, partition.edge_mutations(graph, key, &MutationTrack::Nucleotide)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(AncestralOutputMaps {
    root_sequence,
    node_sequences,
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
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.augur_node_sequence(key))
    })
    .collect();
  let edge_subs = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
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
  progress: &dyn ProgressSink,
) -> Result<AaNodeData, Report> {
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

  progress.check_cancelled()?;
  progress.report("AA ancestral reconstruction", 0.75, "");

  let params = MarginalPartitionParams {
    dense: ancestral_args.dense,
    include_leaves: ancestral_args.include_leaves || ancestral_args.reconstruct_tip_states,
    impute_missing_data: ancestral_args.impute_missing_data || ancestral_args.reconstruct_tip_states,
    sample_from_profile: ancestral_args.sample_from_profile,
    seed: ancestral_args.seed,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  // Reconstruct one CDS partition at a time and consume its result before building the next. A
  // marginal partition holds per-edge probability vectors over the ~20-symbol amino-acid alphabet, so
  // keeping every CDS partition resident at once made peak memory scale with the CDS count. The RNG is
  // created once and passed to each partition so sampled reconstruction draws in a fixed CDS order,
  // independent of how many partitions are resident.
  let mut rng = get_random_number_generator(params.seed);
  let mut aa_node_data = AaNodeData::default();
  for (index, cds) in cdses.iter().enumerate() {
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

    let plan = PartitionPlan {
      name: cds.clone(),
      alphabet: recon_alphabet.clone(),
      gtr_model: aa_model.gtr_model,
      sequences,
      annotation: annotations.get(cds).cloned(),
      reference_override: aa_root_sequences.get(cds).cloned(),
    };

    let reconstructed = reconstruct_marginal_partition(graph, index, plan, &params, names, branch_lengths, &mut rng)?;
    let guard = &reconstructed.partition;

    if let Some(annotation) = &reconstructed.annotation
      && let Some(cds_len) = annotation_cds_nuc_length(annotation)
    {
      let aa_len = i64::try_from(guard.sequence_length())?;
      if 3 * aa_len != cds_len {
        return make_error!(
          "Translated alignment for CDS '{}' has {aa_len} amino acids ({} nucleotides), which does not match \
           the annotated CDS length of {cds_len} nucleotides. Check that the annotation matches the translations.",
          reconstructed.name,
          3 * aa_len
        );
      }
    }

    let cds_data = collect_aa_cds_node_data(
      graph,
      guard,
      &reconstructed.name,
      names,
      reconstructed.reference_override.as_ref(),
    )?;
    aa_node_data.add_cds(&reconstructed.name, cds_data, reconstructed.annotation.clone());

    if let Some(aa_seq_template) = aa_fasta_template {
      write_aa_partition_sequences(graph, guard, names, &reconstructed.name, aa_seq_template)?;
    }
  }

  Ok(aa_node_data)
}

fn write_aa_partition_sequences(
  graph: &Graph,
  partition: &AncestralPartition,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  name: &str,
  template: &str,
) -> Result<(), Report> {
  let path = translation_path(template, name);
  let file = create_file_or_stdout(path)?;
  let mut writer = FastaWriter::new(file);

  for node in graph.get_nodes() {
    let node_key = node.read_arc().key();
    let node_name = names[&node_key]
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);
    let seq = partition.augur_node_sequence(node_key);
    writer.write(&node_name, &None, &seq)?;
  }

  Ok(())
}
