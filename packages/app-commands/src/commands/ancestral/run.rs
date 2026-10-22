use crate::commands::ancestral::aa_node_data::{
  cds_output_paths, read_aa_root_sequences, read_gff3_annotations, selected_cdses, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::{TreetimeAncestralArgs, ancestral_params};
use crate::commands::shared::alignment::pair_alignment;
use crate::commands::shared::ancestral_trees::{
  AncestralOutputMaps, AncestralPartitionTrees, AncestralTrees, write_ancestral_trees,
};
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::reconstructed_fasta::NucFastaWriters;
use crate::commands::shared::sequence_inputs::{
  AlignmentSource, SequenceInputArgs, TreeAndSources, ancestral_input, read_tree_and_sources,
};
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::augur_node_data_ancestral::AncestralNodeSequences;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs, output_unavailable};
use eyre::Report;
use itertools::izip;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::aa::{AaNodeData, AaParams, CdsInput, reconstruct_aa};
use treetime::ancestral::attach::sanitize_to_alphabet;
use treetime::ancestral::params::MethodAncestral;
use treetime::ancestral::pipeline;
use treetime::cancel::Cancel;
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_warn;
use treetime::seq::alignment::{LeafSequences, SequenceSource};
use treetime::seq::gap_fill::{GapFill, apply_gap_fill};
use treetime::seq::sink::{SeqItem, SeqSink, SeqStream, SeqTrack};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_io::fasta::{FastaWriter, fasta_read_file};
use treetime_io::nwk::NwkParse;
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
    args.alignment.alignment.len(),
  )?;

  let sequence_args = SequenceInputArgs {
    alignment: &args.alignment,
    tree: args.tree(),
    tree_dialect: args.tree_dialect.dialect(),
    alphabet_args: &args.alphabet_args,
    gap_fill_args: &args.gap_fill_args,
    ignore_missing_alns: args.missing_sequences.ignore_missing_alns,
  };
  let TreeAndSources { tree, sources } = read_tree_and_sources(&sequence_args, cancel, stages, log)?;
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = tree;
  let (sources, descs): (Vec<SequenceSource>, BTreeMap<usize, _>) = sources
    .into_iter()
    .zip(1..)
    .map(|(AlignmentSource { source, descs, .. }, partition)| (source, (partition, descs)))
    .unzip();
  let input = ancestral_input(graph, branch_lengths, sources);
  let branch_lengths = input.branch_lengths();

  let topology_order =
    args
      .topology_order
      .resolve_topology_order(&input.graph, &names, None, args.tree_dialect.dialect())?;

  let resolved = args.resolve_outputs()?;
  let ambiguity = &args.ambiguity_outputs;
  ambiguity.warn_unwritten(&resolved, log);
  if args.method_anc == MethodAncestral::Parsimony && ambiguity.is_requested() {
    progress_warn!(
      log,
      "--ambiguity-threshold has no effect with --method-anc=parsimony: Fitch parsimony computes no posterior to \
       report ambiguous states from."
    );
  }
  let aa_plan = plan_aa_reconstructions(args, &resolved, log)?;
  let mut seq_sink = AncestralSeqSink::new(&resolved, descs, ambiguity)?;

  let random_step = (args.sample_from_profile != SampleMode::Argmax).then_some("Sampling from the profile");
  let seed = args.seed_args.resolve(random_step, log);
  let params = ancestral_params(args, seed);

  let output = pipeline::run(
    &params,
    input,
    Some(SeqStream::new(&mut seq_sink, &names)),
    cancel,
    stages,
    log,
  )
  .map_err(|err| err.into_report())?;
  let mut graph = output.graph;
  let models = partition_models(&output.results, None);
  let mut results = output.results;
  let AncestralSeqSink {
    fasta,
    mut node_sequences,
    ..
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
  let aa_ambiguity = ambiguity.threshold_for(OutputSelection::ReconstructedAaFasta);
  let aa_result = aa_plan
    .map(|plan| run_aa_reconstructions(args, plan, &aa_inputs, aa_ambiguity, cancel, stages, log))
    .transpose()?;

  let maps = results
    .iter_mut()
    .map(|result| {
      let Some(sequences) = result.sequences.take() else {
        return make_internal_error!("Ancestral reconstruction returned no reconstructed sequences");
      };
      Ok(AncestralOutputMaps {
        root_sequence: sequences.root_sequence,
        edge_mutations: UnknownMutationFilter::new(result.alphabet.unknown(), args.tip_sequences.report_ambiguous)
          .reported_edge_mutations(&graph, sequences.edge_mutations)?,
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  write_gtr_output(&resolved, models, "no GTR model was fitted (use --model=infer)", log)?;

  let trees = AncestralTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    amino_acids: aa_result.as_ref(),
  };
  let partitions = izip!(&results, maps)
    .map(|(result, maps)| AncestralPartitionTrees {
      number: result.number,
      alphabet: &result.alphabet,
      maps,
      node_sequences: node_sequences
        .remove(&result.number)
        .map(|node_sequences| AncestralNodeSequences {
          node_sequences,
          alignment_length: result.length,
          ambiguous_char: result.alphabet.unknown(),
          mask: &result.mask,
        }),
    })
    .collect();
  write_ancestral_trees(&trees, partitions, &resolved, CommandKind::Ancestral, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
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

struct AncestralSeqSink {
  fasta: Option<NucFastaWriters>,
  fasta_ambiguity: Option<AmbiguityThreshold>,
  node_sequences: BTreeMap<usize, BTreeMap<GraphNodeKey, Seq>>,
  node_data_ambiguity: Option<AmbiguityThreshold>,
}

impl AncestralSeqSink {
  fn new(
    resolved: &ResolvedOutputs,
    descs: BTreeMap<usize, BTreeMap<GraphNodeKey, Option<String>>>,
    ambiguity: &AmbiguityOutputs,
  ) -> Result<Self, Report> {
    let fasta_paths = resolved.partition_paths(OutputSelection::ReconstructedNucFasta);
    let fasta = (!fasta_paths.is_empty())
      .then(|| NucFastaWriters::create(fasta_paths, descs))
      .transpose()?;
    let node_sequences: BTreeMap<_, _> = resolved
      .partition_paths(OutputSelection::AugurNodeData)
      .into_keys()
      .map(|partition| (partition, BTreeMap::new()))
      .collect();
    Ok(Self {
      fasta_ambiguity: fasta
        .as_ref()
        .and_then(|_| ambiguity.threshold_for(OutputSelection::ReconstructedNucFasta)),
      fasta,
      node_data_ambiguity: (!node_sequences.is_empty())
        .then(|| ambiguity.threshold_for(OutputSelection::AugurNodeData))
        .flatten(),
      node_sequences,
    })
  }
}

impl SeqSink for AncestralSeqSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    let SeqTrack::Nuc(partition) = item.track else {
      return make_internal_error!("Nucleotide reconstructed FASTA sink received an amino-acid track");
    };
    if let (Some(fasta), true) = (self.fasta.as_mut(), item.emitted) {
      match self.fasta_ambiguity {
        Some(_) => fasta.write(partition, item.key, item.name, &item.seq_with_ambiguity())?,
        None => fasta.write(partition, item.key, item.name, item.seq)?,
      }
    }
    if let Some(node_sequences) = self.node_sequences.get_mut(&partition) {
      let seq = match self.node_data_ambiguity {
        Some(_) => item.seq_with_ambiguity().into_owned(),
        None => item.seq.clone(),
      };
      node_sequences.insert(item.key, seq);
    }
    Ok(())
  }

  fn ambiguity(&self) -> Option<AmbiguityThreshold> {
    self.fasta_ambiguity.or(self.node_data_ambiguity)
  }
}

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  plan: AaPlan<'_>,
  inputs: &AaRunInputs<'_>,
  ambiguity: Option<AmbiguityThreshold>,
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
  let gap_fill_mode = ancestral_args.gap_fill_args.gap_fill();

  let aa_root_sequences = read_aa_root_sequences(ancestral_args.aa_root_sequence.as_deref(), &cdses, &recon_alphabet)?;

  cancel.check()?;
  stages.report("AA ancestral reconstruction", 0.75, "");

  let params = AaParams {
    dense: ancestral_args.representation.dense,
    include_leaves: ancestral_args.include_leaves,
    impute_missing_data: ancestral_args.impute_missing_data,
    sample_from_profile: ancestral_args.sample_from_profile,
    seed: inputs.seed,
    ignore_missing_alns: ancestral_args.missing_sequences.ignore_missing_alns,
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

  let mut seq_sink = fasta_paths.map(|paths| AaFastaSink::new(paths, ambiguity));

  let cds_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry> = cdses
    .iter()
    .filter_map(|cds| annotations.get(cds).map(|entry| (cds.clone(), entry.clone())))
    .collect();

  let node_data = reconstruct_aa(
    inputs.graph,
    inputs.branch_lengths,
    &params,
    cds_inputs,
    seq_sink.as_mut().map(|sink| SeqStream::new(sink, inputs.names)),
    cancel,
    log,
  )?;
  if let Some(sink) = seq_sink {
    sink.finish()?;
  }
  let node_data = UnknownMutationFilter::new(recon_alphabet.unknown(), ancestral_args.tip_sequences.report_ambiguous)
    .reported_aa_node_data(inputs.graph, node_data)?;
  Ok((node_data, cds_annotations))
}

struct AaPlan<'a> {
  translations: &'a str,
  annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
  cdses: Vec<String>,
  fasta_paths: Option<BTreeMap<String, PathBuf>>,
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
  Ok(pair_alignment(sequences, &path, inputs.graph, inputs.names, log).sequences)
}

struct AaRunInputs<'a> {
  graph: &'a Graph,
  names: &'a NodeNames,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  seed: u64,
}

struct AaFastaSink {
  paths: BTreeMap<String, PathBuf>,
  ambiguity: Option<AmbiguityThreshold>,
  pending: Option<(String, BTreeMap<GraphNodeKey, (String, Seq)>)>,
}

impl AaFastaSink {
  fn new(paths: BTreeMap<String, PathBuf>, ambiguity: Option<AmbiguityThreshold>) -> Self {
    Self {
      paths,
      ambiguity,
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
    for (name, seq) in sequences.values() {
      writer.write(name, None, seq)?;
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
    let seq = match self.ambiguity {
      Some(_) => item.seq_with_ambiguity().into_owned(),
      None => item.seq.clone(),
    };
    self
      .pending
      .get_or_insert_with(|| (cds.to_owned(), BTreeMap::new()))
      .1
      .insert(item.key, (item.name.to_owned(), seq));
    Ok(())
  }

  fn ambiguity(&self) -> Option<AmbiguityThreshold> {
    self.ambiguity
  }
}
