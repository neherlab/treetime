use crate::commands::ancestral::aa_node_data::{
  cds_output_paths, read_aa_root_sequences, read_gff3_annotations, selected_cdses, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::{TreetimeAncestralArgs, ancestral_params};
use crate::commands::shared::alignment::pair_alignment;
use crate::commands::shared::ancestral_trees::{AncestralOutputMaps, AncestralTrees, write_ancestral_trees};
use crate::commands::shared::gtr_output::write_gtr_output;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::sequence_inputs::{AncestralReadInputs, SequenceInputArgs, read_nwk_fasta};
use app_output::augur_node_data_ancestral::AncestralNodeSequences;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs, output_unavailable};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::aa::{AaNodeData, AaParams, CdsInput, reconstruct_aa};
use treetime::ancestral::attach::sanitize_to_alphabet;
use treetime::ancestral::pipeline;
use treetime::cancel::Cancel;
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_warn;
use treetime::seq::alignment::LeafSequences;
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime::seq::sink::{SeqItem, SeqSink, SeqTrack};
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::{FastaWriter, fasta_read_file};
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;
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

  let sequence_args = SequenceInputArgs {
    alignment: &args.alignment,
    tree: args.tree(),
    tree_dialect: args.tree_dialect.dialect(),
    alphabet_args: &args.alphabet_args,
    gap_fill_args: &args.gap_fill_args,
    ignore_missing_alns: args.ignore_missing_alns,
  };
  let AncestralReadInputs { input, descs, .. } = read_nwk_fasta(&sequence_args, cancel, stages, log)?;
  let names = input.names();
  let branch_lengths = input.branch_lengths();

  let topology_order =
    args
      .topology_order
      .resolve_topology_order(&input.graph, &names, None, args.tree_dialect.dialect())?;

  let resolved = args.resolve_outputs()?;
  let aa_plan = plan_aa_reconstructions(args, &resolved, log)?;
  let mut seq_sink = AncestralSeqSink::new(&resolved, names.clone(), descs)?;

  let random_step = (args.sample_from_profile != SampleMode::Argmax).then_some("Sampling from the profile");
  let seed = args.seed_args.resolve(random_step, log);
  let params = ancestral_params(args, seed);

  let alphabet = input.alphabet.clone();
  let output =
    pipeline::run(&params, input, Some(&mut seq_sink), cancel, stages, log).map_err(|err| err.into_report())?;
  let mut graph = output.graph;
  let AncestralSeqSink {
    fasta, node_sequences, ..
  } = seq_sink;
  if let Some(fasta) = fasta {
    fasta.finish()?;
  }

  let aa_inputs = AaRunInputs {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    seed,
  };
  let aa_result = aa_plan
    .map(|plan| run_aa_reconstructions(args, plan, &aa_inputs, cancel, stages, log))
    .transpose()?;

  let maps = AncestralOutputMaps {
    root_sequence: output.root_sequence,
    edge_mutations: UnknownMutationFilter::new(output.ambiguous_char, args.report_ambiguous)
      .reported_edge_mutations(&graph, output.edge_mutations)?,
  };

  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  write_gtr_output(
    &resolved,
    output.gtr.as_ref().map(|gtr| (gtr, output.model_name)),
    "no GTR model was fitted (use --model=infer or --gtr-iterations)",
    log,
  )?;

  let trees = AncestralTrees {
    graph: &graph,
    alphabet: &alphabet,
    names: &names,
    branch_lengths: &branch_lengths,
    maps: &maps,
    amino_acids: aa_result.as_ref(),
  };
  let node_sequences = node_sequences.map(|node_sequences| AncestralNodeSequences {
    node_sequences,
    alignment_length: output.sequence_length,
    ambiguous_char: output.ambiguous_char,
    mask: &output.mask,
  });
  write_ancestral_trees(&trees, node_sequences, &resolved, CommandKind::Ancestral, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
}

struct AaRunInputs<'a> {
  graph: &'a Graph,
  names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  seed: u64,
}

fn plan_aa_reconstructions<'a>(
  args: &'a TreetimeAncestralArgs,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<Option<AaPlan<'a>>, Report> {
  let aa_fasta = resolved.non_tree_outputs.get(&OutputSelection::ReconstructedAaFasta);
  let Some(translations) = &args.translations else {
    if let Some(file) = aa_fasta {
      output_unavailable(
        OutputSelection::ReconstructedAaFasta,
        file,
        "no amino-acid translations were given (use --translations)",
        log,
      )?;
    }
    return Ok(None);
  };
  let annotations = read_gff3_annotations(args.annotation.as_deref(), &args.cdses)?;
  let cdses = selected_cdses(&args.cdses, &annotations);
  let fasta_paths = aa_fasta
    .map(|file| planned_aa_fasta_paths(resolved, &file.path, &cdses))
    .transpose()?;
  Ok(Some(AaPlan {
    translations,
    annotations,
    cdses,
    fasta_paths,
  }))
}

fn planned_aa_fasta_paths(
  resolved: &ResolvedOutputs,
  template: &Path,
  cdses: &[String],
) -> Result<BTreeMap<String, PathBuf>, Report> {
  let paths = cds_output_paths(template, cdses)?;
  let flag = OutputSelection::ReconstructedAaFasta.flag_name();
  resolved.ensure_unique_with_expansion(
    OutputSelection::ReconstructedAaFasta,
    paths
      .iter()
      .map(|(cds, path)| (format!("{flag} for CDS '{cds}'"), path.as_path())),
  )?;
  Ok(paths)
}

struct AaPlan<'a> {
  translations: &'a str,
  annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  cdses: Vec<String>,
  fasta_paths: Option<BTreeMap<String, PathBuf>>,
}

struct AncestralSeqSink {
  fasta: Option<FastaWriter>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  descs: BTreeMap<GraphNodeKey, Option<String>>,
  node_sequences: Option<BTreeMap<GraphNodeKey, Seq>>,
}

impl AncestralSeqSink {
  fn new(
    resolved: &ResolvedOutputs,
    names: BTreeMap<GraphNodeKey, Option<String>>,
    descs: BTreeMap<GraphNodeKey, Option<String>>,
  ) -> Result<Self, Report> {
    let fasta = resolved
      .path(OutputSelection::ReconstructedNucFasta)
      .map(FastaWriter::create)
      .transpose()?;
    Ok(Self {
      fasta,
      names,
      descs,
      node_sequences: resolved
        .path(OutputSelection::AugurNodeData)
        .is_some()
        .then(BTreeMap::new),
    })
  }
}

impl SeqSink for AncestralSeqSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    if let (Some(writer), true) = (self.fasta.as_mut(), item.emitted) {
      let name = self.names[&item.key].as_deref();
      let desc = self.descs.get(&item.key).cloned().flatten();
      writer.write(name.unwrap_or(""), desc.as_deref(), item.seq)?;
    }
    if let Some(node_sequences) = self.node_sequences.as_mut() {
      node_sequences.insert(item.key, item.seq.clone());
    }
    Ok(())
  }
}

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  plan: AaPlan<'_>,
  inputs: &AaRunInputs<'_>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>), Report> {
  let AaPlan {
    translations,
    annotations,
    cdses,
    fasta_paths,
  } = plan;
  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let aa_model = ancestral_args.aa_model.resolve();
  let recon_alphabet = Alphabet::new(aa_model.alphabet)?;
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();

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
        sequences: read_cds_translations(
          translations,
          cds,
          &read_alphabet,
          &recon_alphabet,
          gap_fill_mode,
          inputs,
          log,
        )?,
        annotation: annotations.get(cds).cloned(),
        reference_override: aa_root_sequences.get(cds).cloned(),
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut seq_sink = fasta_paths.map(|paths| AaFastaSink::new(paths, inputs.names.clone()));

  let cds_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry> = cdses
    .iter()
    .filter_map(|cds| annotations.get(cds).map(|entry| (cds.clone(), entry.clone())))
    .collect();

  let node_data = reconstruct_aa(
    inputs.graph,
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
  inputs: &AaRunInputs<'_>,
  log: &dyn LogSink,
) -> Result<LeafSequences, Report> {
  let path = translation_path(translations, cds);
  let mut sequences = fasta_read_file(&path, read_alphabet)?;
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
  let paths = [path];
  Ok(pair_alignment(sequences, &paths, inputs.graph, inputs.names, log).sequences)
}

struct AaFastaSink {
  paths: BTreeMap<String, PathBuf>,
  names: BTreeMap<GraphNodeKey, Option<String>>,
  pending: Option<(String, BTreeMap<GraphNodeKey, Seq>)>,
}

impl AaFastaSink {
  fn new(paths: BTreeMap<String, PathBuf>, names: BTreeMap<GraphNodeKey, Option<String>>) -> Self {
    Self {
      paths,
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
    let Some(path) = self.paths.get(&cds) else {
      return make_internal_error!("No amino-acid FASTA path was planned for CDS '{cds}'");
    };
    let mut writer = FastaWriter::create(path)?;
    for (key, seq) in &sequences {
      writer.write(&node_name_or_key(*key, self.names[key].as_deref()), None, seq)?;
    }
    writer.finish()
  }
}

impl SeqSink for AaFastaSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    let SeqTrack::Aa(cds) = item.track else {
      return make_internal_error!("Amino-acid reconstructed FASTA sink received a nucleotide track");
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
