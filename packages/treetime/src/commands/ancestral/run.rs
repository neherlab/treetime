use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::sanitize_to_alphabet;
use crate::ancestral::multi::{MarginalPartitionParams, PartitionPlan, reconstruct_marginal_partition};
use crate::ancestral::pipeline::{self, AncestralInput, AncestralParams, AncestralPartition};
use crate::commands::ancestral::aa_node_data::{
  AaNodeData, annotation_cds_nuc_length, collect_aa_cds_node_data, read_aa_root_sequences, read_gff3_annotations,
  template_has_cds_placeholder, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json_with_aa;
use crate::commands::ancestral::result::{
  AncestralNodeOut, AncestralOutputMaps, AncestralResult, AugurOutputMaps, EdgeOut,
};
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_ancestral_tree_outputs;
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::io::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::{BranchTopology, PartitionBranchOps};
use crate::progress::ProgressSink;
use crate::seq::gap_fill::apply_gap_fill;
use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, mutation_event_strings};
use eyre::Report;
use itertools::Itertools;
use log::{info, warn};
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::nwk::{CommentProviders, NodeCommentProvider};
use treetime_io::nwk::{NwkParse, nwk_read_file};
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use treetime_utils::sync::random::get_random_number_generator;

pub fn run_ancestral_reconstruction(
  ancestral_args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(ancestral_args.alphabet_args.alphabet.unwrap_or_default())?;

  validate_aa_args(
    &ancestral_args.translations,
    &ancestral_args.cdses,
    &ancestral_args.annotation,
    &ancestral_args.aa_root_sequence,
  )?;

  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let mut aln = if ancestral_args.alignment.alignment.is_empty() {
    info!("Reading input fasta from standard input");
    let mut reader = FastaReader::new(open_stdin()?, &alphabet);
    let mut records = Vec::new();
    loop {
      let mut record = FastaRecord::default();
      reader.read(&mut record)?;
      if record.is_empty() {
        break;
      }
      records.push(record);
    }
    records
  } else {
    read_many_fasta(&ancestral_args.alignment.alignment, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  progress.check_cancelled()?;
  progress.report("Parsing tree", 0.1, "");
  let NwkParse {
    graph,
    confidences,
    names,
    branch_lengths: branch_lengths_opt,
  } = nwk_read_file(ancestral_args.tree())?;
  let topology_order = ancestral_args
    .topology_order
    .resolve_topology_order(&graph, &names, None)?;

  let resolved = ancestral_args.resolve_outputs()?;
  let mut output_fasta = if resolved
    .non_tree_outputs
    .contains_key(&OutputSelection::ReconstructedNucFasta)
  {
    let path = &resolved.non_tree_outputs[&OutputSelection::ReconstructedNucFasta];
    Some(FastaWriter::new(create_file_or_stdout(path)?))
  } else {
    None
  };

  let params = AncestralParams {
    method: ancestral_args.method_anc,
    model: ancestral_args.model_args.model,
    dense: ancestral_args.dense,
    include_leaves: ancestral_args.include_leaves || ancestral_args.reconstruct_tip_states,
    impute_missing_data: ancestral_args.impute_missing_data || ancestral_args.reconstruct_tip_states,
    gtr_iterations: ancestral_args.gtr_iterations,
    site_specific_gtr: ancestral_args.site_specific_gtr,
    seed: ancestral_args.seed,
    sample_from_profile: ancestral_args.sample_from_profile,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  // Every reconstruction consumer reads its node label from the `names` map from the parse and its
  // edge branch length from the parsed `branch_lengths` value map.
  // Ancestral never renames or re-lengths after parse, so this map stays accurate at any later point.
  // The pipeline and the amino-acid path derive
  // the `f64` profile map (missing weight resolved to `0.0`) for the marginal passes from it, while
  // the `Option<f64>` map is read directly by the output writers and gather.

  let input = AncestralInput {
    graph,
    alphabet,
    sequences: aln,
  };

  let result = pipeline::run(
    &params,
    input,
    &names,
    &branch_lengths_opt,
    |key, seq| {
      if let Some(ref mut writer) = output_fasta {
        let name = names[&key].as_deref().unwrap_or("");
        // Descriptions originate only on leaf FASTA records and are written during partition init,
        // which runs inside `pipeline::run` below. This writer is the ancestral-reconstruction
        // consumer whose snapshot of node labels is taken before that init; at that point no node
        // carries a description, so every reconstructed record is emitted without one.
        writer.write(name, &None, seq)
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

  let aa_node_data = if let Some(translations) = &ancestral_args.translations {
    Some(run_aa_reconstructions(
      ancestral_args,
      translations,
      aa_fasta_template.as_deref(),
      &result.output.graph,
      &names,
      &branch_lengths_opt,
      progress,
    )?)
  } else {
    // Prerequisite gating: reconstructed AA FASTA needs --translations. An explicit per-file flag
    // is a hard error; the same output reached via selection or `--output-selection=all` is skipped.
    if aa_fasta_template.is_some() {
      if ancestral_args.output_reconstructed_aa_fasta.is_some() {
        return make_error!("--output-reconstructed-aa-fasta requires --translations");
      }
      warn!("Skipping reconstructed amino-acid FASTA output: --translations not provided");
    }
    None
  };

  let pipeline::AncestralOutputFull { output, partition } = result;
  let pipeline::AncestralOutput {
    mut graph,
    gtr,
    model_name,
    mask,
    ..
  } = output;

  // Gather the per-node/per-edge sequence and mutation values off the pipeline-local partition into
  // plain value maps the output writers consume. This is the only place that reads sequences and
  // mutations from the partition; the tree, node-data, and Newick-comment writers read the maps
  // instead. Node and edge keys stay stable through topology ordering, so gathering before it is
  // bit-identical.
  let tree_maps = gather_ancestral_output_maps(&graph, partition.as_ref())?;
  let augur_maps = if resolved.non_tree_outputs.contains_key(&OutputSelection::AugurNodeData) {
    gather_augur_output_maps_opt(&graph, partition.as_ref())?
  } else {
    None
  };

  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;
  progress.report("Writing output", 0.9, "");

  // Gather the per-node name/confidence and per-edge branch length off the ordered tree into keyed
  // value maps the output writers consume. The writers still read sequences and model metadata from
  // the graph data slot; these maps carry the name, input-branch-support, and branch-length values.
  let nodes: BTreeMap<GraphNodeKey, AncestralNodeOut> = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      let confidence = confidences.get(&key).copied().flatten();
      (
        key,
        AncestralNodeOut {
          name: names[&key].clone(),
          confidence,
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      (
        key,
        EdgeOut {
          branch_length: branch_lengths_opt[&key],
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
      write_augur_node_data_json_with_aa(&graph, augur_maps, &mask, &node_names, aa_node_data.as_ref(), path)?;
    }
    info!("Wrote augur node data JSON to {}", path.display());
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, model_name);
        write_gtr_json(&gtr_output, path)?;
      },
      None if ancestral_args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Use --model=infer or --gtr-iterations.");
      },
      None => warn!("Skipping GTR output: no GTR model was fitted (use --model=infer or --gtr-iterations)"),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    write_tree_for_partition(
      &graph,
      &nodes,
      &branch_lengths,
      &tree_maps,
      aa_node_data.as_ref(),
      &resolved,
      partition.as_ref(),
    )?;
  }

  progress.report("Done", 1.0, "");
  Ok(AncestralResult { graph, nodes, edges })
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
  // gathered per-edge mutation map; the partition selects only which provider to use.
  match partition {
    Some(AncestralPartition::Sparse(_) | AncestralPartition::Dense(_)) => {
      let provider = AncestralMutationComments {
        edge_mutations: &maps.edge_mutations,
        graph,
      };
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

/// Gather the per-node nucleotide sequences, root sequence, and per-edge nucleotide mutations the tree
/// writers read off the ancestral partition.
pub(crate) fn gather_ancestral_output_maps(
  graph: &Graph,
  partition: Option<&AncestralPartition>,
) -> Result<AncestralOutputMaps, Report> {
  let Some(partition) = partition else {
    return Ok(AncestralOutputMaps::default());
  };
  match partition {
    AncestralPartition::Fitch(partition) => gather_tree_output_maps(graph, partition),
    AncestralPartition::Sparse(partition) => gather_tree_output_maps(graph, &partition.readout()),
    AncestralPartition::Dense(partition) => gather_tree_output_maps(graph, &partition.readout()),
  }
}

fn gather_tree_output_maps(graph: &Graph, partition: &dyn PartitionBranchOps) -> Result<AncestralOutputMaps, Report> {
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
      Ok((key, partition.edge_mutations(graph, key, MutationTrack::Nucleotide)?))
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
  let maps = match partition {
    AncestralPartition::Fitch(partition) => gather_augur_output_maps(graph, partition)?,
    AncestralPartition::Sparse(partition) => gather_augur_output_maps(graph, partition)?,
    AncestralPartition::Dense(partition) => gather_augur_output_maps(graph, partition)?,
  };
  Ok(Some(maps))
}

/// Gather the augur node-data sequences and substitutions from one partition.
pub(crate) fn gather_augur_output_maps(
  graph: &Graph,
  partition: &dyn AugurNodeDataJsonAncestralPartition,
) -> Result<AugurOutputMaps, Report> {
  let sequence_length = partition.sequence_length();
  let ambiguous_char = partition.ambiguous_char();
  let root_sequence = partition.root_sequence(graph)?;
  let node_sequences = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.node_sequence(key))
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

/// Newick/Nexus node-comment provider that reads the gathered per-edge nucleotide mutation map.
///
/// Mirrors `MutationCommentProvider`, but reads mutations from a value map instead of the partition, so
/// the tree writers no longer touch the partition during serialization.
struct AncestralMutationComments<'a> {
  edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  graph: &'a dyn BranchTopology,
}

impl NodeCommentProvider for AncestralMutationComments<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let Some((_parent_key, edge_key)) = self.graph.node_parent(key)? else {
      return Ok(BTreeMap::new());
    };
    let mut mutations = self.edge_mutations[&edge_key].clone();
    if mutations.is_empty() {
      return Ok(BTreeMap::new());
    }
    mutations.sort_by_key(|mutation| match &mutation.event {
      MutationEvent::Substitution(substitution) => substitution.pos(),
      MutationEvent::Insertion(segment) | MutationEvent::Deletion(segment) => segment.range.0,
    });
    let mutations = mutations
      .iter()
      .map(|mutation| mutation_event_strings(&mutation.event))
      .collect::<Result<Vec<_>, _>>()?
      .into_iter()
      .flatten()
      .join(",");
    Ok(btreemap! {
      "mutations".to_owned() => mutations,
    })
  }
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
    let mut sequences = read_many_fasta(&[&path], &read_alphabet)?;
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
    let guard = reconstructed.partition.as_ref();

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
  partition: &dyn AugurNodeDataJsonAncestralPartition,
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
    let seq = partition.node_sequence(node_key);
    writer.write(&node_name, &None, &seq)?;
  }

  Ok(())
}
