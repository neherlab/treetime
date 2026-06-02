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
use crate::gtr::get_gtr::write_gtr_json;
use crate::make_error;
use crate::partition::traits::MutationCommentProvider;
use crate::payload::ancestral::GraphAncestral;
use crate::progress::ProgressSink;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::{info, warn};
use treetime_graph::node::Named;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::graph::write_graph_files_with_options;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};

pub fn run_ancestral_reconstruction(
  ancestral_args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  let outdir = &ancestral_args.output.outdir;
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

  let output_fasta = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut output_fasta = FastaWriter::new(output_fasta);

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
      let name = node.name.as_deref().unwrap_or("");
      let desc = &node.desc;
      output_fasta.write(name, desc, seq)
    },
    progress,
  )?;

  let aa_node_data = if let Some(translations) = &ancestral_args.translations {
    Some(run_aa_reconstructions(
      ancestral_args,
      translations,
      &result.output.graph,
      progress,
    )?)
  } else {
    None
  };

  progress.report("Writing output", 0.9, "");
  let augur_node_data_path = ancestral_args
    .output_augur_node_data
    .clone()
    .unwrap_or_else(|| outdir.join("ancestral.augur-node-data.json"));

  if let Some(gtr) = &result.output.gtr {
    write_gtr_json(gtr, result.output.model_name, outdir, None)?;
  }

  match &result.partition {
    Some(AncestralPartition::Fitch(partition)) => {
      let graph_options = ancestral_args.output.graph_write_options(&result.output.graph)?;
      let guard = partition.read_arc();
      write_augur_node_data_json_with_aa(
        &result.output.graph,
        &*guard,
        &result.output.mask,
        aa_node_data.as_ref(),
        &augur_node_data_path,
      )?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      write_graph_files_with_options(
        outdir,
        "annotated_tree",
        &result.output.graph,
        &CommentProviders::new(),
        &graph_options,
      )?;
    },
    Some(AncestralPartition::Sparse(partition)) => {
      let graph_options = ancestral_args.output.graph_write_options(&result.output.graph)?;
      let guard = partition.read_arc();
      write_augur_node_data_json_with_aa(
        &result.output.graph,
        &*guard,
        &result.output.mask,
        aa_node_data.as_ref(),
        &augur_node_data_path,
      )?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_graph_files_with_options(
        outdir,
        "annotated_tree",
        &result.output.graph,
        &providers,
        &graph_options,
      )?;
    },
    Some(AncestralPartition::Dense(partition)) => {
      let graph_options = ancestral_args.output.graph_write_options(&result.output.graph)?;
      let guard = partition.read_arc();
      write_augur_node_data_json_with_aa(
        &result.output.graph,
        &*guard,
        &result.output.mask,
        aa_node_data.as_ref(),
        &augur_node_data_path,
      )?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_graph_files_with_options(
        outdir,
        "annotated_tree",
        &result.output.graph,
        &providers,
        &graph_options,
      )?;
    },
    None => {
      let graph_options = ancestral_args.output.graph_write_options(&result.output.graph)?;
      write_graph_files_with_options(
        outdir,
        "annotated_tree",
        &result.output.graph,
        &CommentProviders::new(),
        &graph_options,
      )?;
    },
  }

  progress.report("Done", 1.0, "");
  Ok(AncestralResult {
    graph: result.output.graph,
    gtr: result.output.gtr,
    model_name: result.output.model_name,
  })
}

fn run_aa_reconstructions(
  ancestral_args: &TreetimeAncestralArgs,
  translations: &str,
  graph: &GraphAncestral,
  progress: &dyn ProgressSink,
) -> Result<AaNodeData, Report> {
  // AA translations can contain stop codons (`*`), so read with the stop-inclusive alphabet and then
  // fold out-of-alphabet characters into the unknown state for empirical models that lack a stop.
  let read_alphabet = Alphabet::new(AlphabetName::Aa)?;
  let aa_model = ancestral_args.aa_model.resolve();
  let recon_alphabet = Alphabet::new(aa_model.alphabet)?;
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();

  let annotations = read_gff3_annotations(ancestral_args.annotation.as_deref(), &ancestral_args.cdses)?;

  // When --cdses is omitted, derive the CDS set from the annotation.
  let cdses: Vec<String> = if ancestral_args.cdses.is_empty() {
    annotations.keys().cloned().collect()
  } else {
    ancestral_args.cdses.clone()
  };

  if let Some(aa_seq_template) = &ancestral_args.output_aa_sequences {
    if cdses.len() > 1 && !template_has_cds_placeholder(aa_seq_template) {
      return make_error!(
        "--output-aa-sequences template needs a CDS placeholder when reconstructing multiple CDSes, \
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

  // Write per-CDS reconstructed FASTA when --output-aa-sequences is given.
  if let Some(aa_seq_template) = &ancestral_args.output_aa_sequences {
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
