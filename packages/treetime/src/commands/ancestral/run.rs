use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::pipeline::{self, AncestralInput, AncestralParams, AncestralPartition};
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json;
use crate::commands::ancestral::result::AncestralResult;
use crate::gtr::get_gtr::write_gtr_json;
use crate::partition::traits::MutationCommentProvider;
use crate::progress::ProgressSink;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::graph::write_graph_files_with;
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
      let guard = partition.read_arc();
      write_augur_node_data_json(&result.output.graph, &*guard, &result.output.mask, &augur_node_data_path)?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      write_graph_files_with(outdir, "annotated_tree", &result.output.graph, &CommentProviders::new())?;
    },
    Some(AncestralPartition::Sparse(partition)) => {
      let guard = partition.read_arc();
      write_augur_node_data_json(&result.output.graph, &*guard, &result.output.mask, &augur_node_data_path)?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_graph_files_with(outdir, "annotated_tree", &result.output.graph, &providers)?;
    },
    Some(AncestralPartition::Dense(partition)) => {
      let guard = partition.read_arc();
      write_augur_node_data_json(&result.output.graph, &*guard, &result.output.mask, &augur_node_data_path)?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      let provider = MutationCommentProvider::new(&*guard, &result.output.graph);
      let providers = CommentProviders::new().with(&provider);
      write_graph_files_with(outdir, "annotated_tree", &result.output.graph, &providers)?;
    },
    None => {
      write_graph_files_with(outdir, "annotated_tree", &result.output.graph, &CommentProviders::new())?;
    },
  }

  progress.report("Done", 1.0, "");
  Ok(AncestralResult {
    graph: result.output.graph,
    gtr: result.output.gtr,
    model_name: result.output.model_name,
  })
}
