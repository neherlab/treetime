use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences, create_fitch_partition};
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
use crate::ancestral::mask::create_mask;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::augur_node_data::write_augur_node_data_json;
use crate::commands::ancestral::result::AncestralResult;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr, write_gtr_json};
use crate::gtr::refinement::refine_gtr_iterative;
use crate::make_error;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::fitch::PartitionFitch;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::traits::{HasGtr, MutationCommentProvider};
use crate::payload::ancestral::GraphAncestral;
use crate::progress::ProgressSink;
use crate::seq::alignment::get_common_length;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::graph::write_graph_files_with;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};
use treetime_utils::sync::random::get_random_number_generator;

pub fn run_ancestral_reconstruction(
  ancestral_args: &TreetimeAncestralArgs,
  progress: &dyn ProgressSink,
) -> Result<AncestralResult, Report> {
  let input_fastas = &ancestral_args.alignment.alignment;
  let tree = &ancestral_args.tree;
  let model_name = &ancestral_args.model_args.model;
  let reconstruct_tip_states = &ancestral_args.reconstruct_tip_states;
  let method_anc = &ancestral_args.method_anc;
  let dense = &ancestral_args.dense;
  let gtr_iterations = &ancestral_args.gtr_iterations;
  let outdir = &ancestral_args.output.outdir;
  let site_specific_gtr = &ancestral_args.site_specific_gtr;

  if *site_specific_gtr {
    return make_error!(
      "--site-specific-gtr is not yet integrated into the ancestral reconstruction pipeline. \
       The mathematical core (GTRSiteSpecific) is implemented but partition system wiring is pending."
    );
  }

  if ancestral_args.sample_from_profile != SampleMode::Argmax && *method_anc != MethodAncestral::Marginal {
    return make_error!(
      "--sample-from-profile={:?} requires --method-anc=marginal. Posterior sampling is only defined \
       for marginal reconstruction; {:?} has no posterior profile to sample. Use --method-anc=marginal, \
       or --sample-from-profile=argmax.",
      ancestral_args.sample_from_profile,
      method_anc
    );
  }

  let dense = dense.unwrap_or_else(infer_dense);

  let alphabet = Alphabet::new(ancestral_args.alphabet_args.alphabet.unwrap_or_default())?;

  let gap_fill_mode = ancestral_args.gap_fill_args.effective_gap_fill();

  let mut aln = if input_fastas.is_empty() {
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
    read_many_fasta(input_fastas, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let output_fasta = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut output_fasta = FastaWriter::new(output_fasta);

  progress.check_cancelled()?;
  progress.report("Parsing tree", 0.1, "");
  let graph: GraphAncestral = nwk_read_file(tree)?;

  // Per-position ambiguity mask and augur node data output path are shared by
  // all reconstruction paths below. The mask is derived from the alignment, so
  // it is computed once here while `aln` and `alphabet` are still available.
  let alignment_length = get_common_length(&aln)?;
  let mask = create_mask(&aln, alignment_length, &alphabet);
  let augur_node_data_path = ancestral_args
    .output_augur_node_data
    .clone()
    .unwrap_or_else(|| outdir.join("ancestral.augur-node-data.json"));

  // Posterior-profile sampling for marginal reconstruction. `Argmax` (default) is deterministic and
  // ignores the RNG; `Root`/`All` draw from the per-node posterior. The RNG is seeded from `--seed`
  // for reproducibility, or from entropy when unset.
  let sample_mode = ancestral_args.sample_from_profile;
  let mut rng = get_random_number_generator(ancestral_args.seed);

  match method_anc {
    MethodAncestral::Parsimony => {
      progress.check_cancelled()?;
      progress.report("Fitch parsimony", 0.3, "");
      let partitions_parsimony = vec![Arc::new(RwLock::new(PartitionFitch {
        index: 0,
        alphabet,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }))];

      if !partitions_parsimony.is_empty() {
        compress_sequences(&graph, &partitions_parsimony, &aln)?;

        ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, &partitions_parsimony, |node, seq| {
          let name = node.payload.name.as_deref().unwrap_or("");
          let desc = &node.payload.desc;
          output_fasta.write(name, desc, seq)
        })?;
      }

      progress.report("Writing output", 0.9, "");
      write_augur_node_data_json(
        &graph,
        &*partitions_parsimony[0].read_arc(),
        &mask,
        &augur_node_data_path,
      )?;
      info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
      write_graph_files_with(outdir, "annotated_tree", &graph, &CommentProviders::new())?;

      progress.report("Done", 1.0, "");
      Ok(AncestralResult {
        graph,
        gtr: None,
        model_name: *model_name,
      })
    },
    MethodAncestral::Marginal => {
      progress.check_cancelled()?;
      progress.report("Inferring GTR model", 0.2, "");
      if !dense {
        let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
        let gtr = match *model_name {
          GtrModelName::Infer => infer_gtr_fitch(&fitch, &graph)?,
          _ => get_gtr_by_name(*model_name)?,
        };
        log_gtr(&gtr, *model_name);
        let partition = fitch.into_marginal_sparse(gtr, &graph)?;
        let partitions = vec![Arc::new(RwLock::new(partition))];

        progress.check_cancelled()?;
        progress.report("Marginal reconstruction", 0.4, "");
        update_marginal(&graph, &partitions)?;

        if *gtr_iterations > 0 && *model_name == GtrModelName::Infer {
          refine_gtr_iterative(&graph, &partitions[0], *gtr_iterations, None, 1.0, None, false)?;
        }

        progress.check_cancelled()?;
        progress.report("Reconstructing sequences", 0.6, "");
        ancestral_reconstruction_marginal(
          &graph,
          *reconstruct_tip_states,
          &partitions,
          sample_mode,
          &mut rng,
          |node, seq| {
            let name = node.name.as_deref().unwrap_or("");
            let desc = &node.desc;
            output_fasta.write(name, desc, seq)
          },
        )?;

        progress.report("Writing output", 0.9, "");
        let gtr = partitions[0].read_arc().gtr().clone();
        write_gtr_json(&gtr, *model_name, outdir, None)?;

        let partition_guard = partitions[0].read_arc();
        write_augur_node_data_json(&graph, &*partition_guard, &mask, &augur_node_data_path)?;
        info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
        let provider = MutationCommentProvider::new(&*partition_guard, &graph);
        let providers = CommentProviders::new().with(&provider);
        write_graph_files_with(outdir, "annotated_tree", &graph, &providers)?;

        progress.report("Done", 1.0, "");
        Ok(AncestralResult {
          graph,
          gtr: Some(gtr),
          model_name: *model_name,
        })
      } else if *model_name == GtrModelName::Infer {
        let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
        let gtr = infer_gtr_fitch(&fitch, &graph)?;
        log_gtr(&gtr, *model_name);
        let partition = fitch.into_marginal_dense(gtr);
        let partitions = vec![Arc::new(RwLock::new(partition))];

        progress.check_cancelled()?;
        progress.report("Marginal reconstruction", 0.4, "");
        initialize_marginal(&graph, &partitions, &aln)?;
        update_marginal(&graph, &partitions)?;

        if *gtr_iterations > 0 {
          refine_gtr_iterative(&graph, &partitions[0], *gtr_iterations, None, 1.0, None, false)?;
        }

        progress.check_cancelled()?;
        progress.report("Reconstructing sequences", 0.6, "");
        ancestral_reconstruction_marginal(
          &graph,
          *reconstruct_tip_states,
          &partitions,
          sample_mode,
          &mut rng,
          |node, seq| {
            let name = node.name.as_deref().unwrap_or("");
            let desc = &node.desc;
            output_fasta.write(name, desc, seq)
          },
        )?;

        progress.report("Writing output", 0.9, "");
        let gtr = partitions[0].read_arc().gtr().clone();
        write_gtr_json(&gtr, *model_name, outdir, None)?;

        let partition_guard = partitions[0].read_arc();
        write_augur_node_data_json(&graph, &*partition_guard, &mask, &augur_node_data_path)?;
        info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
        let provider = MutationCommentProvider::new(&*partition_guard, &graph);
        let providers = CommentProviders::new().with(&provider);
        write_graph_files_with(outdir, "annotated_tree", &graph, &providers)?;

        progress.report("Done", 1.0, "");
        Ok(AncestralResult {
          graph,
          gtr: Some(gtr),
          model_name: *model_name,
        })
      } else {
        let length = get_common_length(&aln)?;
        let gtr = get_gtr_by_name(*model_name)?;
        write_gtr_json(&gtr, *model_name, outdir, None)?;
        log_gtr(&gtr, *model_name);
        let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);
        let partitions = vec![Arc::new(RwLock::new(partition))];

        progress.check_cancelled()?;
        progress.report("Marginal reconstruction", 0.4, "");
        initialize_marginal(&graph, &partitions, &aln)?;
        update_marginal(&graph, &partitions)?;

        progress.check_cancelled()?;
        progress.report("Reconstructing sequences", 0.6, "");
        ancestral_reconstruction_marginal(
          &graph,
          *reconstruct_tip_states,
          &partitions,
          sample_mode,
          &mut rng,
          |node, seq| {
            let name = node.name.as_deref().unwrap_or("");
            let desc = &node.desc;
            output_fasta.write(name, desc, seq)
          },
        )?;

        progress.report("Writing output", 0.9, "");
        let gtr = partitions[0].read_arc().gtr().clone();

        let partition_guard = partitions[0].read_arc();
        write_augur_node_data_json(&graph, &*partition_guard, &mask, &augur_node_data_path)?;
        info!("Wrote augur node data JSON to {}", augur_node_data_path.display());
        let provider = MutationCommentProvider::new(&*partition_guard, &graph);
        let providers = CommentProviders::new().with(&provider);
        write_graph_files_with(outdir, "annotated_tree", &graph, &providers)?;

        progress.report("Done", 1.0, "");
        Ok(AncestralResult {
          graph,
          gtr: Some(gtr),
          model_name: *model_name,
        })
      }
    },
  }
}
