use crate::commands::ancestral::aa_node_data::{
  read_aa_root_sequences, read_gff3_annotations, selected_cdses, template_has_cds_placeholder, translation_path,
  validate_aa_args,
};
use crate::commands::ancestral::args::{TreetimeAncestralArgs, ancestral_params};
use crate::commands::shared::alignment::sequence_descriptions;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::EdgeMutationCommentProvider;
use app_output::ancestral_result::{AncestralNodeOut, AncestralOutputMaps, AugurOutputMaps};
use app_output::ancestral_tree_output::write_ancestral_tree_outputs;
use app_output::augur_node_data_ancestral::write_augur_node_data_json_with_aa;
use app_output::output_plan::OutputSelection;
use eyre::Report;
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::aa::{AaNodeData, AaParams, CdsInput, reconstruct_aa};
use treetime::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
use treetime::ancestral::mask::create_mask;
use treetime::ancestral::params::MethodAncestral;
use treetime::ancestral::pipeline;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
use treetime::seq::gap_fill::apply_gap_fill;
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime::{progress_info, progress_warn};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaReader, FastaWriter, read_many_fasta, read_many_fasta_path};
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub fn run_ancestral_reconstruction(
  args: &TreetimeAncestralArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  validate_aa_args(
    args.translations.as_deref(),
    &args.cdses,
    args.annotation.as_deref(),
    args.aa_root_sequence.as_deref(),
  )?;

  let AncestralReadInputs {
    mut input,
    descs,
    confidences,
  } = read_nwk_fasta(args, cancel, stages, log)?;
  let names = input.names();
  let branch_lengths = input.branch_lengths();

  let topology_order = args.topology_order.resolve_topology_order(&input.graph, &names, None)?;

  let resolved = args.resolve_outputs()?;
  let fasta = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedNucFasta)
    .map(|path| Ok::<_, Report>(FastaWriter::new(create_file_or_stdout(path)?)))
    .transpose()?;
  let mut seq_sink = AncestralSeqSink {
    fasta,
    names: names.clone(),
    descs,
    node_sequences: resolved
      .non_tree_outputs
      .contains_key(&OutputSelection::AugurNodeData)
      .then(BTreeMap::new),
  };

  let params = ancestral_params(args);

  let output = pipeline::run(&params, &input, &mut seq_sink, cancel, stages, log).map_err(|err| err.into_report())?;

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
      stages,
      log,
    )?)
  } else {
    if aa_fasta_template.is_some() {
      if args.output_reconstructed_aa_fasta.is_some() {
        return make_error!("--output-reconstructed-aa-fasta requires --translations");
      }
      progress_warn!(
        log,
        "Skipping reconstructed amino-acid FASTA output: --translations not provided"
      );
    }
    None
  };

  let pipeline::AncestralOutput {
    gtr,
    model_name,
    method,
    mask,
    sequence_length,
    ambiguous_char,
    root_sequence,
    edge_mutations,
  } = output;
  let maps = AncestralOutputMaps {
    root_sequence,
    edge_mutations,
  };
  let augur_maps = seq_sink.node_sequences.map(|node_sequences| AugurOutputMaps {
    node_sequences,
    sequence_length,
    ambiguous_char,
  });

  topology_order.apply(&mut input.graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

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

  let aa_node_data = aa_result.as_ref().map(|(node_data, _)| node_data);
  let empty_aa_annotations = BTreeMap::new();
  let aa_annotations = aa_result
    .as_ref()
    .map_or(&empty_aa_annotations, |(_, annotations)| annotations);

  if let (Some(path), Some(augur_maps)) = (
    resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData),
    &augur_maps,
  ) {
    write_augur_node_data_json_with_aa(
      &input.graph,
      &maps,
      augur_maps,
      &mask,
      &names,
      aa_node_data,
      aa_annotations,
      path,
    )?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::builder().gtr(gtr).model_name(model_name).build();
        write_gtr_json(&gtr_output, path)?;
      },
      None if args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Use --model=infer or --gtr-iterations.");
      },
      None => progress_warn!(
        log,
        "Skipping GTR output: no GTR model was fitted (use --model=infer or --gtr-iterations)"
      ),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    write_ancestral_trees(
      &input.graph,
      &nodes,
      &branch_lengths,
      &maps,
      aa_node_data,
      aa_annotations,
      &resolved,
      method,
    )?;
  }

  stages.report("Done", 1.0, "");
  Ok(())
}

struct AncestralSeqSink {
  fasta: Option<FastaWriter>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  descs: BTreeMap<String, Option<String>>,
  node_sequences: Option<BTreeMap<GraphNodeKey, Seq>>,
}

impl SeqSink for AncestralSeqSink {
  fn on_topology(&mut self, _graph: &Graph) -> Result<(), Report> {
    Ok(())
  }

  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    if let (Some(writer), true) = (self.fasta.as_mut(), item.emitted) {
      let name = self.names[&item.key].as_deref();
      let desc = name.and_then(|name| self.descs.get(name)).cloned().flatten();
      writer.write(name.unwrap_or(""), &desc, item.seq)?;
    }
    if let Some(node_sequences) = self.node_sequences.as_mut() {
      node_sequences.insert(item.key, item.seq.clone());
    }
    Ok(())
  }
}

fn read_nwk_fasta(
  args: &TreetimeAncestralArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralReadInputs, Report> {
  let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let mut aln = if args.alignment.alignment.is_empty() {
    progress_info!(log, "Reading input fasta from standard input");
    let reader = FastaReader::new(open_stdin()?, &alphabet);
    read_many_fasta(reader)?
  } else {
    read_many_fasta_path(&args.alignment.alignment, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  let descs = sequence_descriptions(&aln);

  cancel.check()?;
  stages.report("Parsing tree", 0.1, "");
  let parse = nwk_read_file(args.tree())?;

  let confidences = parse.confidences();

  let names = parse.names();
  let aln = aln.into_iter().map(AlignmentRecord::from).collect();
  let aln = complete_alignment_for_leaves(&parse.graph, aln, &alphabet, args.ignore_missing_alns, &names, log)?;
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);

  let graph = parse.graph;
  let nodes = node_seq_inputs(&graph, &names, aln);
  let edges = parse
    .branch_lengths
    .into_iter()
    .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
    .collect();
  let input = AncestralInput {
    graph,
    nodes,
    edges,
    alphabet,
    mask,
  };
  Ok(AncestralReadInputs {
    input,
    descs,
    confidences,
  })
}

struct AncestralReadInputs {
  input: AncestralInput,
  descs: BTreeMap<String, Option<String>>,
  confidences: BTreeMap<GraphNodeKey, Option<f64>>,
}

fn write_ancestral_trees(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
  aa_annotations: &BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  resolved: &app_output::output_plan::ResolvedOutputs,
  method: MethodAncestral,
) -> Result<(), Report> {
  let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, graph);
  let providers = match method {
    MethodAncestral::Marginal => CommentProviders::new().with(&provider),
    MethodAncestral::Parsimony | MethodAncestral::Joint => CommentProviders::new(),
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

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  translations: &str,
  aa_fasta_template: Option<&str>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>), Report> {
  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let aa_model = ancestral_args.aa_model.resolve();
  let recon_alphabet = Alphabet::new(aa_model.alphabet)?;
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();

  let annotations = read_gff3_annotations(ancestral_args.annotation.as_deref(), &ancestral_args.cdses)?;

  let cdses = selected_cdses(&ancestral_args.cdses, &annotations);

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
  stages.report("AA ancestral reconstruction", 0.75, "");

  let params = AaParams {
    dense: ancestral_args.dense,
    include_leaves: ancestral_args.include_leaves,
    impute_missing_data: ancestral_args.impute_missing_data,
    sample_from_profile: ancestral_args.sample_from_profile,
    seed: ancestral_args.seed,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  let mut cds_inputs = Vec::with_capacity(cdses.len());
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
      progress_warn!(
        log,
        "CDS '{cds}': mapped {sanitized} out-of-alphabet amino-acid characters (e.g. stop '*') to '{}'.",
        char::from(recon_alphabet.unknown())
      );
    }

    cds_inputs.push(CdsInput {
      name: cds.clone(),
      alphabet: recon_alphabet.clone(),
      gtr_model: aa_model.gtr_model,
      sequences: sequences.into_iter().map(AlignmentRecord::from).collect(),
      annotation: annotations.get(cds).cloned(),
      reference_override: aa_root_sequences.get(cds).cloned(),
    });
  }

  let mut seq_sink = aa_fasta_template.map(|template| AaFastaSink::new(template.to_owned(), names.clone()));

  let cds_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry> = cdses
    .iter()
    .filter_map(|cds| annotations.get(cds).map(|entry| (cds.clone(), entry.clone())))
    .collect();

  let node_data = reconstruct_aa(
    graph,
    names,
    branch_lengths,
    &params,
    cds_inputs,
    seq_sink.as_mut().map(|sink| -> &mut dyn SeqSink { sink }),
    cancel,
    log,
  )?;
  if let Some(sink) = seq_sink {
    sink.finish()?;
  }
  Ok((node_data, cds_annotations))
}

struct AaFastaSink {
  template: String,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  pending: Option<(String, BTreeMap<GraphNodeKey, Seq>)>,
}

impl AaFastaSink {
  fn new(template: String, names: BTreeMap<GraphNodeKey, Option<String>>) -> Self {
    Self {
      template,
      names,
      pending: None,
    }
  }

  fn finish(mut self) -> Result<(), Report> {
    self.write_pending()
  }

  fn write_pending(&mut self) -> Result<(), Report> {
    let Some((cds, sequences)) = self.pending.take() else {
      return Ok(());
    };
    let mut writer = FastaWriter::new(create_file_or_stdout(translation_path(&self.template, &cds))?);
    for (key, seq) in &sequences {
      let name = self.names[key]
        .as_deref()
        .map_or_else(|| format!("node_{}", key.0), str::to_owned);
      writer.write(&name, &None, seq)?;
    }
    Ok(())
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
    if self.pending.as_ref().is_some_and(|(pending, _)| pending != cds) {
      self.write_pending()?;
    }
    self
      .pending
      .get_or_insert_with(|| (cds.to_owned(), BTreeMap::new()))
      .1
      .insert(item.key, item.seq.clone());
    Ok(())
  }
}
