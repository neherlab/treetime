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
use app_output::gtr::write_gtr_json;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{OutputSelection, ResolvedOutputs};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::aa::{AaNodeData, AaParams, CdsInput, reconstruct_aa};
use treetime::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
use treetime::ancestral::mask::create_mask;
use treetime::ancestral::params::MethodAncestral;
use treetime::ancestral::pipeline;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput};
use treetime::gtr::gtr::GTR;
use treetime::make_error;
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::{LogSink, StageSink};
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
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
  let mut seq_sink = AncestralSeqSink::new(&resolved, names.clone(), descs)?;

  let random_step = (args.sample_from_profile != SampleMode::Argmax).then_some("Sampling from the profile");
  let seed = args.seed_args.resolve(random_step, log);
  let params = ancestral_params(args, seed);

  let output = pipeline::run(&params, &input, &mut seq_sink, cancel, stages, log).map_err(|err| err.into_report())?;
  let AncestralSeqSink {
    fasta, node_sequences, ..
  } = seq_sink;
  if let Some(fasta) = fasta {
    fasta.finish()?;
  }

  let aa_inputs = AaRunInputs {
    graph: &input.graph,
    names: &names,
    branch_lengths: &branch_lengths,
    seed,
  };
  let aa_result = optional_aa_reconstructions(args, &resolved, &aa_inputs, cancel, stages, log)?;

  let maps = AncestralOutputMaps {
    root_sequence: output.root_sequence,
    edge_mutations: UnknownMutationFilter::new(output.ambiguous_char, args.report_ambiguous)
      .reported_edge_mutations(&input.graph, output.edge_mutations)?,
  };
  let augur_maps = node_sequences.map(|node_sequences| AugurOutputMaps {
    node_sequences,
    sequence_length: output.sequence_length,
    ambiguous_char: output.ambiguous_char,
  });

  topology_order.apply(&mut input.graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  let nodes = ancestral_node_outputs(&input, &confidences);

  let aa_node_data = aa_result.as_ref().map(|(node_data, _)| node_data);
  let empty_aa_annotations = BTreeMap::new();
  let aa_annotations = aa_result
    .as_ref()
    .map_or(&empty_aa_annotations, |(_, annotations)| annotations);

  if let (Some(path), Some(augur_maps)) = (
    resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData),
    augur_maps,
  ) {
    write_augur_node_data_json_with_aa(
      &input.graph,
      &maps,
      augur_maps,
      &output.mask,
      &names,
      aa_node_data,
      aa_annotations,
      path,
    )?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }

  write_ancestral_gtr(args, &resolved, output.gtr.as_ref(), output.model_name, log)?;

  write_ancestral_trees(
    &input.graph,
    &nodes,
    &branch_lengths,
    &maps,
    aa_node_data,
    aa_annotations,
    &resolved,
    output.method,
    log,
  )?;

  stages.report("Done", 1.0, "");
  Ok(())
}

struct AaRunInputs<'a> {
  graph: &'a Graph,
  names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  seed: u64,
}

fn optional_aa_reconstructions(
  args: &TreetimeAncestralArgs,
  resolved: &ResolvedOutputs,
  inputs: &AaRunInputs<'_>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<Option<(AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>)>, Report> {
  let aa_fasta_template: Option<String> = resolved
    .non_tree_outputs
    .get(&OutputSelection::ReconstructedAaFasta)
    .map(|path| path.to_string_lossy().into_owned());
  let Some(translations) = &args.translations else {
    if aa_fasta_template.is_some() {
      if args.output_reconstructed_aa_fasta.is_some() {
        return make_error!("--output-reconstructed-aa-fasta requires --translations");
      }
      progress_warn!(
        log,
        "Skipping reconstructed amino-acid FASTA output: --translations not provided"
      );
    }
    return Ok(None);
  };
  run_aa_reconstructions(
    args,
    translations,
    aa_fasta_template.as_deref(),
    inputs,
    cancel,
    stages,
    log,
  )
  .map(Some)
}

fn ancestral_node_outputs(
  input: &AncestralInput,
  confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> BTreeMap<GraphNodeKey, AncestralNodeOut> {
  input
    .nodes
    .iter()
    .map(|(key, node)| {
      let out = AncestralNodeOut {
        name: node.name.clone(),
        branch_support: confidences.get(key).copied().flatten(),
      };
      (*key, out)
    })
    .collect()
}

fn write_ancestral_gtr(
  args: &TreetimeAncestralArgs,
  resolved: &ResolvedOutputs,
  gtr: Option<&GTR>,
  model_name: GtrModelName,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) else {
    return Ok(());
  };
  match gtr {
    Some(gtr) => {
      let gtr_output = GtrOutput::builder().gtr(gtr).model_name(model_name).build();
      write_gtr_json(&gtr_output, path)
    },
    None if args.output_gtr.is_some() => {
      make_error!("GTR output requested but no GTR model was fitted. Use --model=infer or --gtr-iterations.")
    },
    None => {
      progress_warn!(
        log,
        "Skipping GTR output: no GTR model was fitted (use --model=infer or --gtr-iterations)"
      );
      Ok(())
    },
  }
}

struct AncestralSeqSink {
  fasta: Option<FastaWriter>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  descs: BTreeMap<String, Option<String>>,
  node_sequences: Option<BTreeMap<GraphNodeKey, Seq>>,
}

impl AncestralSeqSink {
  fn new(
    resolved: &ResolvedOutputs,
    names: BTreeMap<GraphNodeKey, Option<String>>,
    descs: BTreeMap<String, Option<String>>,
  ) -> Result<Self, Report> {
    let fasta = resolved
      .non_tree_outputs
      .get(&OutputSelection::ReconstructedNucFasta)
      .map(|path| Ok::<_, Report>(FastaWriter::new(create_file_or_stdout(path)?)))
      .transpose()?;
    Ok(Self {
      fasta,
      names,
      descs,
      node_sequences: resolved
        .non_tree_outputs
        .contains_key(&OutputSelection::AugurNodeData)
        .then(BTreeMap::new),
    })
  }
}

impl SeqSink for AncestralSeqSink {
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
  resolved: &ResolvedOutputs,
  method: MethodAncestral,
  log: &dyn LogSink,
) -> Result<(), Report> {
  if resolved.tree_outputs.is_empty() {
    return Ok(());
  }
  let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, graph);
  let providers = match method {
    MethodAncestral::Marginal => CommentProviders::new().with(&provider),
    MethodAncestral::Parsimony => CommentProviders::new(),
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
    log,
  )
}

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  translations: &str,
  aa_fasta_template: Option<&str>,
  inputs: &AaRunInputs<'_>,
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
    seed: inputs.seed,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  let cds_inputs = cdses
    .iter()
    .map(|cds| {
      Ok(CdsInput {
        name: cds.clone(),
        alphabet: recon_alphabet.clone(),
        gtr_model: aa_model.gtr_model,
        sequences: read_cds_translations(translations, cds, &read_alphabet, &recon_alphabet, gap_fill_mode, log)?,
        annotation: annotations.get(cds).cloned(),
        reference_override: aa_root_sequences.get(cds).cloned(),
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut seq_sink = aa_fasta_template.map(|template| AaFastaSink::new(template.to_owned(), inputs.names.clone()));

  let cds_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry> = cdses
    .iter()
    .filter_map(|cds| annotations.get(cds).map(|entry| (cds.clone(), entry.clone())))
    .collect();

  let node_data = reconstruct_aa(
    inputs.graph,
    inputs.names,
    inputs.branch_lengths,
    &params,
    cds_inputs,
    seq_sink.as_mut().map(|sink| -> &mut dyn SeqSink { sink }),
    cancel,
    log,
  )?;
  if let Some(sink) = seq_sink {
    sink.finish()?;
  }
  let node_data = UnknownMutationFilter::new(recon_alphabet.unknown(), ancestral_args.report_ambiguous)
    .reported_aa_node_data(inputs.graph, node_data)?;
  Ok((node_data, cds_annotations))
}

fn read_cds_translations(
  translations: &str,
  cds: &str,
  read_alphabet: &Alphabet,
  recon_alphabet: &Alphabet,
  gap_fill_mode: GapFill,
  log: &dyn LogSink,
) -> Result<Vec<AlignmentRecord>, Report> {
  let path = translation_path(translations, cds);
  let mut sequences = read_many_fasta_path(&[&path], read_alphabet)?;
  let mut sanitized = 0_usize;
  for record in &mut sequences {
    let (seq, changed) = sanitize_to_alphabet(&record.seq, recon_alphabet);
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
  Ok(sequences.into_iter().map(AlignmentRecord::from).collect())
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
    writer.finish()
  }
}

impl SeqSink for AaFastaSink {
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
