use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::pipeline::{self, AncestralInput, AncestralParams, AncestralPartition};
use crate::commands::ancestral::aa_node_data::{
  AaNodeData, collect_aa_gene_node_data, read_aa_root_sequences, read_gff3_annotations, translation_path,
  validate_aa_args,
};
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json_with_aa;
use crate::commands::ancestral::result::AncestralResult;
use crate::gtr::get_gtr::{GtrModelName, write_gtr_json};
use crate::partition::traits::MutationCommentProvider;
use crate::progress::ProgressSink;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::graph::write_graph_files_with_options;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::{nwk_read_file, nwk_read_str};
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use treetime_utils::io::fs::read_file_to_string;

pub fn run_ancestral_reconstruction(
  ancestral_args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  let outdir = &ancestral_args.output.outdir;
  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();
  let alphabet = Alphabet::new(ancestral_args.alphabet_args.alphabet.unwrap_or_default())?;

  validate_aa_args(
    &ancestral_args.translations,
    &ancestral_args.genes,
    &ancestral_args.annotation_gff,
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
    Some(run_aa_reconstructions(ancestral_args, translations, progress)?)
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
  progress: &dyn ProgressSink,
) -> Result<AaNodeData, Report> {
  let aa_alphabet = Alphabet::new(AlphabetName::AaNoStop)?;
  let mut aa_node_data = AaNodeData::default();
  let annotations = read_gff3_annotations(ancestral_args.annotation_gff.as_deref(), &ancestral_args.genes)?;
  let aa_root_sequences = read_aa_root_sequences(ancestral_args.aa_root_sequence.as_deref(), &ancestral_args.genes)?;
  let tree = read_file_to_string(&ancestral_args.tree)?;

  for (index, gene) in ancestral_args.genes.iter().enumerate() {
    progress.check_cancelled()?;
    progress.report(
      "AA ancestral reconstruction",
      0.75,
      &format!("{} ({}/{})", gene, index + 1, ancestral_args.genes.len()),
    );

    let path = translation_path(translations, gene);
    let mut sequences = read_many_fasta(&[&path], &aa_alphabet)?;
    let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();
    for record in &mut sequences {
      apply_gap_fill(&mut record.seq, gap_fill_mode, aa_alphabet.gap(), aa_alphabet.unknown());
    }

    let graph = nwk_read_str(&tree)?;
    let input = AncestralInput {
      graph,
      alphabet: aa_alphabet.clone(),
      sequences,
    };
    let params = AncestralParams {
      method: ancestral_args.method_anc,
      model: GtrModelName::Jtt92,
      dense: ancestral_args.dense,
      reconstruct_tip_states: ancestral_args.reconstruct_tip_states,
      gtr_iterations: 0,
      site_specific_gtr: false,
      seed: ancestral_args.seed,
      sample_from_profile: ancestral_args.sample_from_profile,
    };

    let result = pipeline::run(&params, input, |_node, _seq| Ok(()), progress)?;
    let reference_override = aa_root_sequences.get(gene);
    let gene_data = match &result.partition {
      Some(AncestralPartition::Fitch(partition)) => {
        let guard = partition.read_arc();
        collect_aa_gene_node_data(&result.output.graph, &*guard, gene, reference_override)?
      },
      Some(AncestralPartition::Sparse(partition)) => {
        let guard = partition.read_arc();
        collect_aa_gene_node_data(&result.output.graph, &*guard, gene, reference_override)?
      },
      Some(AncestralPartition::Dense(partition)) => {
        let guard = partition.read_arc();
        collect_aa_gene_node_data(&result.output.graph, &*guard, gene, reference_override)?
      },
      None => continue,
    };
    aa_node_data.add_gene(gene, gene_data, annotations.get(gene).cloned());
  }

  Ok(aa_node_data)
}
