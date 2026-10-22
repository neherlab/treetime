use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::shared::gtr_output::{trait_model, write_gtr_output};
use crate::commands::shared::input_warnings::warn_duplicate_names;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::tree_input::read_input_tree;
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeTraits};
use app_output::augur_node_data_traits::write_augur_node_data_traits;
use app_output::command_outputs::CommandKind;
use app_output::output_plan::{OutputSelection, ResolvedOutputs};
use app_output::trait_tables::{write_trait_confidence_csv, write_traits_csv};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::cancel::Cancel;
use treetime::mugration::pipeline::{self, MugrationInput, MugrationOutput, MugrationParams};
use treetime::mugration::trait_pairing::{TraitPairing, pair_trait_rows};
use treetime::progress::RunWarningKind;
use treetime::progress::{LogSink, StageSink};
use treetime::{progress_info, progress_warn};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_io::discrete_states_csv::discrete_attrs_read_file;
use treetime_utils::make_error;

const UNMATCHED_NAMES_SHOWN: usize = 10;

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");
  let resolved = mugration_args.resolve_outputs()?;
  let ambiguity = &mugration_args.ambiguity_outputs;
  ambiguity.warn_unwritten(&resolved, log);
  let parse = read_input_tree(&mugration_args.tree, mugration_args.newick_dialect.dialect(), log)?;
  let names = parse.names;
  let graph: Graph = parse.graph;
  let branch_lengths = parse.branch_lengths;

  let attribute_column = Some(mugration_args.attribute());

  let (rows, _attr_name) = discrete_attrs_read_file::<String>(
    mugration_args.metadata(),
    &mugration_args.metadata_id.metadata_delimiters,
    &mugration_args.metadata_id.metadata_id_columns,
    None,
    attribute_column,
    |s| Ok(s.to_owned()),
  )?;
  let TraitPairing {
    traits,
    observed_values,
    ..
  } = pair_traits(rows, &graph, &names, mugration_args.metadata(), log);

  let weights = if let Some(weights_filepath) = &mugration_args.weights {
    let (map, _) = discrete_attrs_read_file::<f64>(
      weights_filepath,
      &mugration_args.metadata_id.metadata_delimiters,
      &[],
      attribute_column,
      Some("weight"),
      |s| Ok(s.parse::<f64>()?),
    )?;
    let weights = map.into_iter().collect::<BTreeMap<String, f64>>();
    validate_weights(&weights, weights_filepath)?;
    Some(weights)
  } else {
    None
  };

  cancel.check()?;
  stages.report("Mugration inference", 0.3, "");
  let params = MugrationParams {
    missing_data: mugration_args.missing_data.clone(),
    pc: mugration_args.pc,
    missing_weights_threshold: mugration_args.missing_weights_threshold,
    iterations: mugration_args.model_iterations,
    sampling_bias_correction: mugration_args.sampling_bias_correction,
    smooth_initial_pi: mugration_args.smooth_initial_pi,
    filter_uninformative_root: mugration_args.filter_uninformative_root,
  };
  let input = MugrationInput {
    graph,
    traits,
    traits_source: format!("metadata '{}'", mugration_args.metadata().display()),
    observed_values,
    weights,
    branch_lengths: branch_lengths.clone(),
  };
  let mut output = pipeline::run(&params, input, &names, cancel, log).map_err(|err| err.into_report())?;

  let topology_order = mugration_args.topology_order.resolve_topology_order(
    &output.graph,
    &names,
    None,
    mugration_args.newick_dialect.dialect(),
  )?;
  topology_order.apply(&mut output.graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.8, "");

  write_gtr_output(
    &resolved,
    vec![trait_model(
      &output.gtr,
      mugration_args.attribute(),
      output.states.iter().map(ToOwned::to_owned).collect(),
    )],
    "no trait model was fitted",
    log,
  )?;

  write_mugration_trees(
    &output,
    &names,
    &branch_lengths,
    mugration_args.attribute(),
    &resolved,
    ambiguity,
    log,
  )?;

  stages.report("Done", 1.0, "");
  Ok(())
}

pub(crate) fn validate_weights(weights: &BTreeMap<String, f64>, weights_filepath: &Path) -> Result<(), Report> {
  match weights
    .iter()
    .find(|(_, weight)| !(weight.is_finite() && **weight > 0.0))
  {
    Some((state, weight)) => make_error!(
      "Weights file '{}': weight of state '{state}' is {weight}, expected a finite positive number",
      weights_filepath.display()
    ),
    None => Ok(()),
  }
}

pub(crate) fn pair_traits(
  rows: Vec<(String, String)>,
  graph: &Graph,
  names: &NodeNames,
  path: &Path,
  log: &dyn LogSink,
) -> TraitPairing {
  let pairing = pair_trait_rows(rows, graph, names);
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateMetadataNames,
    &format!("The metadata '{}' has more than one row named", path.display()),
    "TreeTime uses the first row of each name.",
    &pairing.duplicate_names,
  );
  if !pairing.unmatched.is_empty() {
    let shown = pairing
      .unmatched
      .iter()
      .take(UNMATCHED_NAMES_SHOWN)
      .map(|(name, _)| name)
      .join(", ");
    let suffix = if pairing.unmatched.len() > UNMATCHED_NAMES_SHOWN {
      "..."
    } else {
      ""
    };
    progress_warn!(
      log,
      "Mugration: {} metadata names not present in tree: {shown}{suffix}",
      pairing.unmatched.len()
    );
  }
  pairing
}

fn write_mugration_trees(
  output: &MugrationOutput,
  names: &NodeNames,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  attribute: &str,
  resolved: &ResolvedOutputs,
  ambiguity: &AmbiguityOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph: &output.graph,
    names,
    divergence_branch_lengths: branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: None,
    dates: None,
    traits: Some(TreeTraits {
      attribute,
      states: &output.states,
      posteriors: &output.posteriors,
      ambiguity: Some(ambiguity),
    }),
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  let Some(tree) = tree_view_for_outputs(&annotated, resolved)? else {
    return Ok(());
  };
  write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Mugration, log)?;
  if let Some(path) = resolved.path(OutputSelection::Traits) {
    write_traits_csv(&tree, path)?;
  }
  if let Some(path) = resolved.path(OutputSelection::TraitProbabilities) {
    write_trait_confidence_csv(&tree, path)?;
  }
  if let Some(path) = resolved.path(OutputSelection::NodeData) {
    write_augur_node_data_traits(&tree, &output.gtr, path)?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }
  Ok(())
}
