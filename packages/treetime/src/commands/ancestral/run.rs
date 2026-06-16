use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::sanitize_to_alphabet;
use crate::ancestral::multi::{MarginalPartitionParams, PartitionPlan, reconstruct_marginal_partitions};
use crate::ancestral::pipeline::{self, AncestralInput, AncestralParams, AncestralPartition};
use crate::commands::ancestral::aa_node_data::{
  AaNodeData, annotation_cds_nuc_length, collect_aa_cds_node_data, read_aa_root_sequences, read_gff3_annotations,
  template_has_cds_placeholder, translation_path, validate_aa_args,
};
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json_with_aa;
use crate::commands::ancestral::result::AncestralResult;
use crate::commands::shared::output::{CommandKind, OutputSelection};
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::traits::MutationCommentProvider;
use crate::payload::ancestral::GraphAncestral;
use crate::progress::ProgressSink;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::{info, warn};
use treetime_graph::node::Named;
use treetime_graph::topology_order::TopologyOrderSpec;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::graph::write_tree_outputs;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};

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
  let graph = nwk_read_file(&ancestral_args.tree)?;
  let topology_order = ancestral_args.topology_order.resolve_topology_order(&graph, None)?;

  let selection: Vec<OutputSelection> = ancestral_args
    .output_selection
    .iter()
    .copied()
    .map(OutputSelection::from)
    .collect();
  let resolved = ancestral_args.output.resolve(
    CommandKind::Ancestral,
    &selection,
    &[
      (
        OutputSelection::AugurNodeData,
        ancestral_args.output_augur_node_data.as_deref(),
      ),
      (OutputSelection::Gtr, ancestral_args.output_gtr.as_deref()),
      (
        OutputSelection::ReconstructedNucFasta,
        ancestral_args.output_reconstructed_nuc_fasta.as_deref(),
      ),
      (
        OutputSelection::ReconstructedAaFasta,
        ancestral_args
          .output_reconstructed_aa_fasta
          .as_deref()
          .map(std::path::Path::new),
      ),
    ],
  )?;

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
    reconstruct_tip_states: ancestral_args.reconstruct_tip_states,
    gtr_iterations: ancestral_args.gtr_iterations,
    site_specific_gtr: ancestral_args.site_specific_gtr,
    seed: ancestral_args.seed,
    sample_from_profile: ancestral_args.sample_from_profile,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  let input = AncestralInput {
    graph,
    alphabet,
    sequences: aln,
  };

  let result = pipeline::run(
    &params,
    input,
    |node, seq| {
      if let Some(ref mut writer) = output_fasta {
        let name = node.name.as_deref().unwrap_or("");
        let desc = &node.desc;
        writer.write(name, desc, seq)
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

  progress.report("Writing output", 0.9, "");

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    match &result.partition {
      Some(AncestralPartition::Fitch(partition)) => {
        let guard = partition.read_arc();
        write_augur_node_data_json_with_aa(
          &result.output.graph,
          &*guard,
          &result.output.mask,
          aa_node_data.as_ref(),
          path,
        )?;
      },
      Some(AncestralPartition::Sparse(partition)) => {
        let guard = partition.read_arc();
        write_augur_node_data_json_with_aa(
          &result.output.graph,
          &*guard,
          &result.output.mask,
          aa_node_data.as_ref(),
          path,
        )?;
      },
      Some(AncestralPartition::Dense(partition)) => {
        let guard = partition.read_arc();
        write_augur_node_data_json_with_aa(
          &result.output.graph,
          &*guard,
          &result.output.mask,
          aa_node_data.as_ref(),
          path,
        )?;
      },
      None => {},
    }
    info!("Wrote augur node data JSON to {}", path.display());
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match result.output.gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, result.output.model_name);
        write_gtr_json(&gtr_output, path)?;
      },
      None if ancestral_args.output_gtr.is_some() => {
        return make_error!("GTR output requested but no GTR model was fitted. Use --model=infer or --gtr-iterations.");
      },
      None => warn!("Skipping GTR output: no GTR model was fitted (use --model=infer or --gtr-iterations)"),
    }
  }

  if !resolved.tree_outputs.is_empty() {
    write_tree_for_partition(&result, &resolved, &topology_order)?;
  }

  progress.report("Done", 1.0, "");
  Ok(AncestralResult {
    graph: result.output.graph,
    gtr: result.output.gtr,
    model_name: result.output.model_name,
  })
}

fn write_tree_for_partition(
  result: &pipeline::AncestralOutputFull,
  resolved: &crate::commands::shared::output::ResolvedOutputs,
  topology_order: &TopologyOrderSpec,
) -> Result<(), Report> {
  let plan = topology_order.plan(&result.output.graph)?;
  let ordered = plan.ordered_graph(&result.output.graph)?;

  match &result.partition {
    Some(AncestralPartition::Sparse(partition)) => {
      let guard = partition.read_arc();
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, None)?;
    },
    Some(AncestralPartition::Dense(partition)) => {
      let guard = partition.read_arc();
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, None)?;
    },
    Some(AncestralPartition::Fitch(_)) | None => {
      write_tree_outputs(&ordered, &resolved.tree_outputs, &CommentProviders::new(), None)?;
    },
  }
  Ok(())
}

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  translations: &str,
  aa_fasta_template: Option<&str>,
  graph: &GraphAncestral,
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

  let mut plans = Vec::with_capacity(cdses.len());
  for cds in &cdses {
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

    plans.push(PartitionPlan {
      name: cds.clone(),
      alphabet: recon_alphabet.clone(),
      gtr_model: aa_model.gtr_model,
      sequences,
      annotation: annotations.get(cds).cloned(),
      reference_override: aa_root_sequences.get(cds).cloned(),
    });
  }

  let params = MarginalPartitionParams {
    dense: ancestral_args.dense,
    reconstruct_tip_states: ancestral_args.reconstruct_tip_states,
    sample_from_profile: ancestral_args.sample_from_profile,
    seed: ancestral_args.seed,
    ignore_missing_alns: ancestral_args.ignore_missing_alns,
  };

  let reconstructed = reconstruct_marginal_partitions(graph, plans, &params)?;

  let mut aa_node_data = AaNodeData::default();
  for partition in &reconstructed {
    let guard = partition.partition.read_arc();

    if let Some(annotation) = &partition.annotation
      && let Some(cds_len) = annotation_cds_nuc_length(annotation)
    {
      let aa_len = i64::try_from(guard.sequence_length())?;
      if 3 * aa_len != cds_len {
        return make_error!(
          "Translated alignment for CDS '{}' has {aa_len} amino acids ({} nucleotides), which does not match \
           the annotated CDS length of {cds_len} nucleotides. Check that the annotation matches the translations.",
          partition.name,
          3 * aa_len
        );
      }
    }

    let cds_data = collect_aa_cds_node_data(graph, &*guard, &partition.name, partition.reference_override.as_ref())?;
    aa_node_data.add_cds(&partition.name, cds_data, partition.annotation.clone());
  }

  if let Some(aa_seq_template) = aa_fasta_template {
    write_aa_sequences(graph, &reconstructed, aa_seq_template)?;
  }

  Ok(aa_node_data)
}

fn write_aa_sequences(
  graph: &GraphAncestral,
  reconstructed: &[crate::ancestral::multi::ReconstructedPartition],
  template: &str,
) -> Result<(), Report> {
  for partition in reconstructed {
    let path = translation_path(template, &partition.name);
    let file = create_file_or_stdout(path)?;
    let mut writer = FastaWriter::new(file);
    let guard = partition.partition.read_arc();

    for node in graph.get_nodes() {
      let node_guard = node.read_arc();
      let payload = node_guard.payload().read_arc();
      let name = payload
        .name()
        .map_or_else(|| format!("node_{}", node_guard.key().0), |n| n.as_ref().to_owned());
      let seq = guard.node_sequence(node_guard.key());
      writer.write(&name, &None, &seq)?;
    }
  }

  Ok(())
}
