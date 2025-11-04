use crate::commands::clock::clock_model::ClockModel;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::seq::div::{OnlyLeaves, calculate_divs};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  keep_root: bool,
  clock_model: &ClockModel,
) -> Result<(), Report> {
  log::info!("# Running timetree inference");

  log::info!("## Calculating divergence distances");
  let divs = calculate_divs(graph, OnlyLeaves(false));
  for node_ref in graph.get_nodes() {
    let mut node = node_ref.write_arc().payload().write_arc();
    if let Some(name) = &node.name {
      if let Some(&div) = divs.get(name) {
        node.div = div;
      }
    }
  }

  log::info!("## Using clock model");
  log::info!("**Clock rate:** {:.6e}", clock_model.clock_rate());

  if !partitions.is_empty() {
    log::info!("## Computing branch distributions from partitions");
    compute_branch_distributions(graph, partitions, clock_model)?;
  } else {
    log::info!("## Creating branch distributions from input lengths");
    create_branch_distributions_from_input_lengths(graph, clock_model)?;
  }

  log::info!("## Propagating distributions backward");
  propagate_distributions_backward(graph)?;

  log::info!("## Propagating distributions forward");
  propagate_distributions_forward(graph)?;

  log::info!("# Timetree inference completed");
  Ok(())
}

fn compute_branch_distributions(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  // In input branch mode, partitions exist but have no edge data
  // Skip branch distribution computation and use input branch lengths directly
  if partitions.iter().any(|p| p.read_arc().get_sequence_length().is_none()) {
    log::debug!("Skipping branch distribution computation: partitions not initialized");
    return Ok(());
  }

  let one_mutation = calculate_one_mutation(partitions);
  let total_sites: usize = partitions
    .iter()
    .map(|p| p.read_arc().get_sequence_length().unwrap_or(0))
    .sum();

  log::info!(
    "Computing branch distributions from {} partition(s) with {} total sites",
    partitions.len(),
    total_sites
  );
  log::debug!("One mutation = {one_mutation:.6e} substitutions/site");

  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.branch_length.unwrap_or(one_mutation);

    log::debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}");

    let contributions = collect_contributions(partitions, edge_key)?;
    let distribution = compute_branch_length_distribution(
      &contributions,
      branch_length,
      one_mutation,
      BRANCH_GRID_SIZE,
      clock_rate,
    )?;

    if let Some(likely_time) = distribution.likely_time() {
      log::debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
    }

    edge.branch_length_distribution = Some(distribution);
  }
  Ok(())
}

fn calculate_one_mutation(partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>]) -> f64 {
  let total_length: usize = partitions
    .iter()
    .map(|part| part.read_arc().get_sequence_length().unwrap_or(0))
    .sum();
  1.0 / total_length as f64
}

fn collect_contributions(
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  edge_key: GraphEdgeKey,
) -> Result<Vec<OptimizationContribution>, Report> {
  partitions
    .iter()
    .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
    .collect()
}

fn create_branch_distributions_from_input_lengths(
  graph: &GraphAncestral,
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length {
      // Convert branch length (substitutions/site) to time duration (years)
      let time_duration = branch_length / clock_rate;
      let distribution = Distribution::point(time_duration, 1.0);
      edge.branch_length_distribution = Some(Arc::new(distribution));
    }
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::commands::ancestral::marginal_unified::run_marginal;
  use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
  use crate::commands::clock::clock_set::ClockSet;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::data::date_constraints::load_date_constraints;
  use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::dates_csv::{DateOrRange, DatesMap};
  use crate::io::fasta::read_many_fasta;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::representation::partition_marginal_dense::PartitionMarginalDense;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
  use crate::seq::div::{OnlyLeaves, calculate_divs};
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_io::json::{JsonPretty, json_write_file};

  #[test]
  fn test_timetree_flu_h3n2_poisson() -> Result<(), Report> {
    const CLOCK_RATE_MU: f64 = 0.0028;
    const SEQUENCE_LENGTH_L: usize = 1400;

    fn dump_graph(graph: &GraphAncestral, filename: &str) -> Result<(), Report> {
      let OUTPUT_DIR = Path::new("../../tmp/test_timetree_flu_h3n2_poisson");
      json_write_file(OUTPUT_DIR.join(filename), graph, JsonPretty(true))?;
      Ok(())
    }

    // Rerooted H3N2 flu tree obtained from:
    //
    // treetime clock \
    //    --tree data/flu/h3n2/20/tree.nwk \
    //    --dates data/flu/h3n2/20/metadata.tsv \
    //    --sequence-length 1400 \
    //    --outdir tmp/clock_test
    //
    // Output: tmp/clock_test/rerooted.newick
    let input_nwk = r#"((A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409:0.00143,A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409:0.00070)NODE_00000161.00:0.01778,(A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416:0.00647,(A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409:0.00210,(A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423:0.00578,(A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409:0.00350,(A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412:0.00428,(A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409:0.00214,((A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409:0.00214,A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416:0.00357)NODE_00000090.98:0.00426,(A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428:0.00761,(A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416:0.00517,((A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409:0.00142,A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409:0.00214)NODE_00000050.93:0.00214,(A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409:0.00285,((A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409:0.00500,A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409:0.00071)NODE_00000010.47:0.00055,(A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409:0.00357,A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409:0.00142)NODE_00000020.45:0.00055)NODE_0000000:0.00071)NODE_00000030.81:0.00286)NODE_00000040.96:0.00486)NODE_00000060.98:0.00184)NODE_00000070.79:0.00199)NODE_00000080.81:0.00074)NODE_00000100.81:0.00071)NODE_00000110.80:0.00055)NODE_00000120.87:0.01374)NODE_00000131.00:0.01010)NODE_00000141.00:0.00066)NODE_00000150.73:0.00408)NODE_0000017:0.00100;"#;

    // From data/flu/h3n2/20/metadata.tsv
    let input_dates: DatesMap = btreemap! {
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.40520192)),
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.83778234)),
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => Some(DateOrRange::YearFraction(2009.48186174)),
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => Some(DateOrRange::YearFraction(2009.5229295)),
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => Some(DateOrRange::YearFraction(2000.13415469)),
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => Some(DateOrRange::YearFraction(2000.68172485)),
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => Some(DateOrRange::YearFraction(2011.98015058)),
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.95550992)),
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.8569473)),
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.84052019)),
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => Some(DateOrRange::YearFraction(2007.48733744)),
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => Some(DateOrRange::YearFraction(2008.15058179)),
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => Some(DateOrRange::YearFraction(2008.86516085)),
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.25735797)),
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.98562628)),
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => Some(DateOrRange::YearFraction(2011.65160849)),
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.11225188)),
    };

    // Python TreeTime v0 baseline results obtained from:
    //
    // treetime \
    //   --tree tmp/clock_test/rerooted.newick \
    //   --keep-root \
    //   --clock-rate 0.0028 \
    //   --dates data/flu/h3n2/20/metadata.tsv \
    //   --branch-length-mode input \
    //   --sequence-len 1400 \
    //   --keep-polytomies \
    //   --outdir tmp/python_treetime_baseline
    //
    // Output: tmp/python_treetime_baseline/dates.tsv
    //
    // The --branch-length-mode input flag makes Python TreeTime use Poisson approximation
    // for branch length distributions (same as this test), enabling direct comparison.
    let expected = btreemap! {
      o!("NODE_0000017") => 1996.895553,
      o!("NODE_00000161.00") => 2002.842738,
      o!("NODE_00000150.73") => 1998.478787,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_00000141.00") => 1998.750301,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_00000131.00") => 2001.455532,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_00000120.87") => 2006.561853,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_00000110.80") => 2006.897741,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_00000100.81") => 2007.200137,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_00000080.81") => 2007.459045,
      o!("NODE_00000090.98") => 2008.623974,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_00000070.79") => 2008.475534,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_00000060.98") => 2009.169225,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_00000040.96") => 2010.619950,
      o!("NODE_00000050.93") => 2011.351158,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_00000030.81") => 2011.528039,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.869954,
      o!("NODE_00000010.47") => 2012.061416,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_00000020.45") => 2012.139155,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
    };

    let graph: GraphAncestral = nwk_read_str(input_nwk)?;

    load_date_constraints(&input_dates, &graph)?;

    create_poisson_branch_distributions(&graph, CLOCK_RATE_MU, SEQUENCE_LENGTH_L, BRANCH_GRID_SIZE)?;
    dump_graph(&graph, "001_after_create_poisson_branch_distributions.json")?;

    propagate_distributions_backward(&graph)?;
    dump_graph(&graph, "002_after_propagate_distributions_backward.json")?;

    propagate_distributions_forward(&graph)?;
    dump_graph(&graph, "003_after_propagate_distributions_forward.json")?;

    let actual: BTreeMap<String, f64> = graph
      .get_nodes()
      .into_iter()
      .filter_map(|node_ref| {
        let node = node_ref.read_arc();
        let payload = node.payload().read_arc();
        let name = payload.name.as_ref()?.clone();
        let time = payload.time?;
        Some((name, time))
      })
      .collect();

    let actual = round_times(actual);
    let expected = round_times(expected);

    // {
    //   let max_diff_actual = expected
    //     .iter()
    //     .filter_map(|(name, &expected_time)| actual.get(name).map(|&actual_time| (expected_time - actual_time).abs()))
    //     .max_by_key(|&x| OrderedFloat(x))
    //     .unwrap_or(0.0);

    //   let max_diff_expected = 0.419;
    //   assert_abs_diff_eq!(max_diff_actual, max_diff_expected, epsilon = 1e-6);
    // }

    // {
    //   let expected_keys: BTreeSet<_> = expected.keys().collect();
    //   let actual_keys: BTreeSet<_> = actual.keys().collect();

    //   let missing_nodes: BTreeSet<_> = expected_keys.difference(&actual_keys).collect();
    //   let extra_nodes: BTreeSet<_> = actual_keys.difference(&expected_keys).collect();

    //   assert_eq!(missing_nodes, btreeset! {}, "Some nodes are missing");
    //   assert_eq!(extra_nodes, btreeset! {}, "Extra nodes found");
    // }

    {
      assert_eq!(&expected.into_iter().collect_vec(), &actual.into_iter().collect_vec());
    }

    Ok(())
  }

  #[test]
  fn test_timetree_flu_h3n2_marginal_sparse() -> Result<(), Report> {
    const SEQUENCE_LENGTH_L: usize = 1407;
    const CLOCK_RATE: f64 = 0.0028;

    fn dump_graph(graph: &GraphAncestral, filename: &str) -> Result<(), Report> {
      let OUTPUT_DIR = Path::new("../../tmp/test_timetree_flu_h3n2_marginal_sparse");
      json_write_file(OUTPUT_DIR.join(filename), graph, JsonPretty(true))?;
      Ok(())
    }

    // Rerooted H3N2 flu tree obtained from:
    //
    // cargo run --bin=treetime -- clock \
    //   --tree data/flu/h3n2/20/tree.nwk \
    //   --dates data/flu/h3n2/20/metadata.tsv \
    //   --outdir tmp/clock_test
    //
    // Output: tmp/clock_test/rerooted.nwk
    let input_nwk = r#"((A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416:0.00647,(A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409:0.0021,(A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423:0.00578,(A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409:0.0035,(A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412:0.00428,(A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409:0.00214,((A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409:0.00214,A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416:0.00357)NODE_0000006:0.00426,(A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428:0.00761,(A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416:0.00517,((A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409:0.00142,A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409:0.00214)NODE_0000014:0.00214,(A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409:0.00285,((A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409:0.005,A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409:0.00071)NODE_0000015:0.00055,(A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409:0.00357,A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409:0.00142)NODE_0000016:0.00055)NODE_0000000:0.00071)NODE_0000001:0.00286)NODE_0000002:0.00486)NODE_0000003:0.00184)NODE_0000004:0.00199)NODE_0000005:0.00074)NODE_0000007:0.00071)NODE_0000008:0.00055)NODE_0000009:0.0137)NODE_0000010:0.0101)NODE_0000011:0.00066)NODE_0000012:0.00393,((A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409:0.00143,A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409:0.0007)NODE_0000013:0.0175):0.000437);"#;

    // From data/flu/h3n2/20/metadata.tsv
    let input_dates: DatesMap = btreemap! {
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.40520192)),
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.83778234)),
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => Some(DateOrRange::YearFraction(2009.48186174)),
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => Some(DateOrRange::YearFraction(2009.5229295)),
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => Some(DateOrRange::YearFraction(2000.13415469)),
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => Some(DateOrRange::YearFraction(2000.68172485)),
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => Some(DateOrRange::YearFraction(2011.98015058)),
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.95550992)),
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.8569473)),
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.84052019)),
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => Some(DateOrRange::YearFraction(2007.48733744)),
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => Some(DateOrRange::YearFraction(2008.15058179)),
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => Some(DateOrRange::YearFraction(2008.86516085)),
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.25735797)),
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.98562628)),
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => Some(DateOrRange::YearFraction(2011.65160849)),
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.11225188)),
    };

    // Python TreeTime v0 baseline results obtained from:
    //
    // treetime \
    //   --tree tmp/clock_test/rerooted.nwk \
    //   --keep-root \
    //   --clock-rate 0.0028 \
    //   --dates data/flu/h3n2/20/metadata.tsv \
    //   --aln data/flu/h3n2/20/aln.fasta \
    //   --sequence-len 1407 \
    //   --keep-polytomies \
    //   --max-iter 0 \
    //   --branch-length-mode marginal \
    //   --time-marginal true \
    //   --outdir tmp/python_treetime_marginal
    //
    // Output: tmp/python_treetime_marginal/dates.tsv
    let expected = btreemap! {
      o!("NODE_0000017") => 1996.090159,
      o!("NODE_0000018") => 1996.988929,
      o!("NODE_0000013") => 2002.843607,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
      o!("NODE_0000012") => 1998.298776,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_0000011") => 1998.764184,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_0000010") => 2001.452891,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_0000009") => 2006.413708,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_0000008") => 2006.685003,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_0000007") => 2007.093618,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_0000005") => 2007.431995,
      o!("NODE_0000006") => 2008.637703,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_0000004") => 2008.400630,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_0000003") => 2009.138143,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_0000002") => 2010.519680,
      o!("NODE_0000014") => 2011.352916,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_0000001") => 2011.456075,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.888447,
      o!("NODE_0000015") => 2012.078684,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_0000016") => 2012.127243,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
    };

    let mut graph: GraphAncestral = nwk_read_str(input_nwk)?;

    load_date_constraints(&input_dates, &graph)?;

    let aln_path = Path::new("../../data/flu/h3n2/20/aln.fasta.xz");
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let aln = read_many_fasta(&[aln_path.to_path_buf()], &alphabet)?;

    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: SEQUENCE_LENGTH_L,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &aln)?;
    dump_graph(&graph, "001_after_compress_sequences.json")?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll>>> = vec![sparse_partition];

    // Run marginal reconstruction
    run_marginal(&graph, &partitions, Some(&aln))?;
    dump_graph(&graph, "002_after_run_marginal.json")?;

    let divs = calculate_divs(&graph, OnlyLeaves(false));
    for node_ref in graph.get_nodes() {
      let mut node = node_ref.write_arc().payload().write_arc();
      if let Some(name) = &node.name {
        if let Some(&div) = divs.get(name) {
          node.div = div;
        }
      }
    }
    dump_graph(&graph, "002_after_calculate_divs.json")?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockOptions::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
    )?;
    dump_graph(&graph, "003_after_clock_model.json")?;

    {
      for node_ref in graph.get_nodes() {
        let mut node = node_ref.write_arc().payload().write_arc();
        if let Some(dist_arc) = &node.time_distribution {
          if let Some(time) = dist_arc.likely_time() {
            node.clock_total = ClockSet::leaf_contribution(Some(time));
          }
        }
      }
    }
    dump_graph(&graph, "004_after_initialize_node_times.json")?;

    run_timetree(&mut graph, &partitions, true, &clock_model)?;
    dump_graph(&graph, "005_after_run_timetree.json")?;

    let actual: BTreeMap<String, f64> = graph
      .get_nodes()
      .into_iter()
      .filter_map(|node_ref| {
        let node = node_ref.read_arc();
        let payload = node.payload().read_arc();
        let name = payload.name.as_ref()?.clone();
        let time = payload.time?;
        Some((name, time))
      })
      .collect();

    let actual = round_times(actual);
    let expected = round_times(expected);

    {
      assert_eq!(&expected.into_iter().collect_vec(), &actual.into_iter().collect_vec());
    }

    Ok(())
  }

  #[test]
  fn test_timetree_flu_h3n2_marginal_dense() -> Result<(), Report> {
    const SEQUENCE_LENGTH_L: usize = 1407;
    const CLOCK_RATE: f64 = 0.0028;

    fn dump_graph(graph: &GraphAncestral, filename: &str) -> Result<(), Report> {
      let OUTPUT_DIR = Path::new("../../tmp/test_timetree_flu_h3n2_marginal_dense");
      json_write_file(OUTPUT_DIR.join(filename), graph, JsonPretty(true))?;
      Ok(())
    }

    // Rerooted H3N2 flu tree obtained from:
    //
    // cargo run --bin=treetime -- clock \
    //   --tree data/flu/h3n2/20/tree.nwk \
    //   --dates data/flu/h3n2/20/metadata.tsv \
    //   --outdir tmp/clock_test
    //
    // Output: tmp/clock_test/rerooted.nwk
    let input_nwk = r#"((A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416:0.00647,(A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409:0.0021,(A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423:0.00578,(A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409:0.0035,(A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412:0.00428,(A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409:0.00214,((A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409:0.00214,A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416:0.00357)NODE_0000006:0.00426,(A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428:0.00761,(A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416:0.00517,((A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409:0.00142,A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409:0.00214)NODE_0000014:0.00214,(A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409:0.00285,((A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409:0.005,A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409:0.00071)NODE_0000015:0.00055,(A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409:0.00357,A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409:0.00142)NODE_0000016:0.00055)NODE_0000000:0.00071)NODE_0000001:0.00286)NODE_0000002:0.00486)NODE_0000003:0.00184)NODE_0000004:0.00199)NODE_0000005:0.00074)NODE_0000007:0.00071)NODE_0000008:0.00055)NODE_0000009:0.0137)NODE_0000010:0.0101)NODE_0000011:0.00066)NODE_0000012:0.00393,((A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409:0.00143,A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409:0.0007)NODE_0000013:0.0175):0.000437);"#;

    // From data/flu/h3n2/20/metadata.tsv
    let input_dates: DatesMap = btreemap! {
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.40520192)),
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.83778234)),
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => Some(DateOrRange::YearFraction(2009.48186174)),
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => Some(DateOrRange::YearFraction(2009.5229295)),
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => Some(DateOrRange::YearFraction(2000.13415469)),
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => Some(DateOrRange::YearFraction(2000.68172485)),
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => Some(DateOrRange::YearFraction(2011.98015058)),
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.95550992)),
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.8569473)),
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.84052019)),
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => Some(DateOrRange::YearFraction(2007.48733744)),
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => Some(DateOrRange::YearFraction(2008.15058179)),
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => Some(DateOrRange::YearFraction(2008.86516085)),
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2012.25735797)),
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => Some(DateOrRange::YearFraction(2011.98562628)),
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => Some(DateOrRange::YearFraction(2003.00273785)),
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => Some(DateOrRange::YearFraction(2011.65160849)),
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => Some(DateOrRange::YearFraction(2013.11225188)),
    };

    // Python TreeTime v0 baseline results obtained from:
    //
    // treetime \
    //   --tree tmp/clock_test/rerooted.nwk \
    //   --keep-root \
    //   --clock-rate 0.0028 \
    //   --dates data/flu/h3n2/20/metadata.tsv \
    //   --aln data/flu/h3n2/20/aln.fasta \
    //   --sequence-len 1407 \
    //   --keep-polytomies \
    //   --max-iter 0 \
    //   --branch-length-mode marginal \
    //   --time-marginal true \
    //   --outdir tmp/python_treetime_marginal
    //
    // Output: tmp/python_treetime_marginal/dates.tsv
    let expected = btreemap! {
      o!("NODE_0000017") => 1996.090159,
      o!("NODE_0000018") => 1996.988929,
      o!("NODE_0000013") => 2002.843607,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
      o!("NODE_0000012") => 1998.298776,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_0000011") => 1998.764184,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_0000010") => 2001.452891,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_0000009") => 2006.413708,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_0000008") => 2006.685003,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_0000007") => 2007.093618,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_0000005") => 2007.431995,
      o!("NODE_0000006") => 2008.637703,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_0000004") => 2008.400630,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_0000003") => 2009.138143,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_0000002") => 2010.519680,
      o!("NODE_0000014") => 2011.352916,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_0000001") => 2011.456075,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.888447,
      o!("NODE_0000015") => 2012.078684,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_0000016") => 2012.127243,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
    };

    let mut graph: GraphAncestral = nwk_read_str(input_nwk)?;

    load_date_constraints(&input_dates, &graph)?;

    let aln_path = Path::new("../../data/flu/h3n2/20/aln.fasta.xz");
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let aln = read_many_fasta(&[aln_path.to_path_buf()], &alphabet)?;

    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: SEQUENCE_LENGTH_L,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll>>> = vec![dense_partition];

    // Run marginal reconstruction using dense representation
    run_marginal(&graph, &partitions, Some(&aln))?;
    dump_graph(&graph, "001_after_run_marginal.json")?;

    let divs = calculate_divs(&graph, OnlyLeaves(false));
    for node_ref in graph.get_nodes() {
      let mut node = node_ref.write_arc().payload().write_arc();
      if let Some(name) = &node.name {
        if let Some(&div) = divs.get(name) {
          node.div = div;
        }
      }
    }
    dump_graph(&graph, "002_after_calculate_divs.json")?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockOptions::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
    )?;
    dump_graph(&graph, "003_after_clock_model.json")?;

    {
      for node_ref in graph.get_nodes() {
        let mut node = node_ref.write_arc().payload().write_arc();
        if let Some(dist_arc) = &node.time_distribution {
          if let Some(time) = dist_arc.likely_time() {
            node.clock_total = ClockSet::leaf_contribution(Some(time));
          }
        }
      }
    }
    dump_graph(&graph, "004_after_initialize_node_times.json")?;

    run_timetree(&mut graph, &partitions, true, &clock_model)?;
    dump_graph(&graph, "005_after_run_timetree.json")?;

    let actual: BTreeMap<String, f64> = graph
      .get_nodes()
      .into_iter()
      .filter_map(|node_ref| {
        let node = node_ref.read_arc();
        let payload = node.payload().read_arc();
        let name = payload.name.as_ref()?.clone();
        let time = payload.time?;
        Some((name, time))
      })
      .collect();

    let actual = round_times(actual);
    let expected = round_times(expected);

    {
      assert_eq!(&expected.into_iter().collect_vec(), &actual.into_iter().collect_vec());
    }

    Ok(())
  }

  fn round_to(x: f64, n: u32) -> f64 {
    let p = 10_f64.powi(n as i32);
    (x * p).round() / p
  }

  fn round_times(times: BTreeMap<String, f64>) -> BTreeMap<String, f64> {
    times
      .into_iter()
      .map(|(name, time)| (name, round_to(time, 3)))
      .collect()
  }

  fn create_poisson_branch_distributions(
    graph: &GraphAncestral,
    mu: f64,
    seq_len: usize,
    n_points: usize,
  ) -> Result<(), Report> {
    let seq_len_f64 = seq_len as f64;

    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();

      if let Some(branch_length) = edge.branch_length {
        let expected_time = branch_length / mu;
        let max_time = 3.0 * expected_time.max(1.0);
        let dt_step = max_time / (n_points - 1) as f64;

        let x = Array1::from_iter((0..n_points).map(|i| i as f64 * dt_step));

        let log_probs = x.mapv(|dt| {
          if dt < 1e-10 {
            f64::NEG_INFINITY
          } else {
            -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln()
          }
        });

        let log_p_max = log_probs
          .iter()
          .copied()
          .filter(|&x| x.is_finite())
          .map(OrderedFloat)
          .max()
          .map_or(0.0, |x| x.0);

        let y = log_probs.mapv(|log_p| {
          if log_p.is_finite() {
            (log_p - log_p_max).exp()
          } else {
            0.0
          }
        });

        let distribution = Distribution::function(x, y)?;
        edge.branch_length_distribution = Some(Arc::new(distribution));
      }
    }

    Ok(())
  }
}
