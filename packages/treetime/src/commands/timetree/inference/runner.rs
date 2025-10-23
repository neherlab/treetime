use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

const BRANCH_GRID_SIZE: usize = 200;

pub fn run_timetree(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  keep_root: bool,
) -> Result<(), Report> {
  log::info!("# Running timetree inference");

  log::info!("## Estimating clock rate from root-to-tip regression");
  let clock_model = estimate_clock_model_with_reroot(
    graph,
    &ClockOptions::default(),
    keep_root,
    &BranchPointOptimizationParams::default(),
  )?;
  log::info!("**Estimated clock rate:** {:.6e}", clock_model.clock_rate());
  log::debug!("Clock rate: {:.6e}", clock_model.clock_rate());

  if !partitions.is_empty() {
    log::info!("## Computing branch distributions from partitions");
    compute_branch_distributions(graph, partitions, &clock_model)?;
  } else {
    log::info!("## Creating branch distributions from input lengths");
    create_branch_distributions_from_input_lengths(graph, &clock_model)?;
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
    return Ok(());
  }

  let one_mutation = calculate_one_mutation(partitions);
  let clock_rate = clock_model.clock_rate();

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let branch_length = edge.branch_length.unwrap_or(one_mutation);

    let contributions = collect_contributions(partitions, edge_key)?;
    let distribution = compute_branch_length_distribution(
      &contributions,
      branch_length,
      one_mutation,
      BRANCH_GRID_SIZE,
      clock_rate,
    )?;
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
  use crate::commands::timetree::data::date_constraints::load_date_constraints;
  use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::io::dates_csv::{DateOrRange, DatesMap};
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use maplit::btreemap;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_timetree_flu_h3n2_poisson() -> Result<(), Report> {
    const CLOCK_RATE_MU: f64 = 0.0028;
    const SEQUENCE_LENGTH_L: usize = 1400;
    const GRID_SIZE: usize = 200;

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
    //   --branch-length-mode input \
    //   --sequence-len 1400 \
    //   --outdir tmp/python_treetime_baseline
    //
    // Output: tmp/python_treetime_baseline/dates.tsv
    //
    // The --branch-length-mode input flag makes Python TreeTime use Poisson approximation
    // for branch length distributions (same as this test), enabling direct comparison.
    let expected = btreemap! {
      o!("NODE_0000017") => 1996.974064,
      o!("NODE_0000012") => 1998.499705,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_0000011") => 1998.763998,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_0000010") => 2001.464619,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_0000009") => 2006.560498,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_0000008") => 2006.896740,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_0000007") => 2007.199495,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_0000005") => 2007.458496,
      o!("NODE_0000006") => 2008.623669,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_0000004") => 2008.475043,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_0000003") => 2009.168874,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_0000002") => 2010.619876,
      o!("NODE_0000014") => 2011.351132,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_0000001") => 2011.528015,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.869945,
      o!("NODE_0000015") => 2012.061411,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_0000016") => 2012.139149,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
      o!("NODE_0000013") => 2002.842738,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
    };

    let graph: GraphAncestral = nwk_read_str(input_nwk)?;

    load_date_constraints(&input_dates, &graph)?;

    create_poisson_branch_distributions(&graph, CLOCK_RATE_MU, SEQUENCE_LENGTH_L, GRID_SIZE)?;

    propagate_distributions_backward(&graph)?;

    propagate_distributions_forward(&graph)?;

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

    assert_eq!(
      json_write_str(&expected, JsonPretty(true))?,
      json_write_str(&actual, JsonPretty(true))?,
    );

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
