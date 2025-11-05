use crate::commands::clock::clock_model::ClockModel;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::utils::{initialize_clock_totals_from_time_distributions, initialize_node_divergences};
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
  clock_model: &ClockModel,
) -> Result<(), Report> {
  log::info!("# Running timetree inference");

  log::info!("## Calculating divergence distances");
  initialize_node_divergences(graph);

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
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::commands::ancestral::marginal_unified::run_marginal;
  use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::data::date_constraints::load_date_constraints;
  use crate::commands::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::commands::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::dates_csv::{DatesMap, read_dates};
  use crate::io::fasta::{FastaRecord, read_many_fasta};
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::representation::partition_marginal_dense::PartitionMarginalDense;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::Path;
  use std::sync::LazyLock;
  use treetime_io::json::{JsonPretty, json_write_file};

  // seqkit fx2tab -l ./data/flu/h3n2/20/aln.fasta.xz | cut -f4 | uniq
  const SEQUENCE_LENGTH: usize = 1407;

  // Output of clock command (see below)
  const CLOCK_RATE: f64 = 0.002826;

  // Rerooted H3N2 flu tree obtained from Treetime v0 (Python), such that we can skip the rerooting step:
  //
  // ./dev/docker/python treetime clock \
  //    --tree data/flu/h3n2/20/tree.nwk \
  //    --dates data/flu/h3n2/20/metadata.tsv \
  //    --sequence-length 1407 \
  //    --outdir tmp/clock_test
  //
  // Output: tmp/clock_test/rerooted.newick
  const INPUT_REROOTED_TREE_NWK: &str = r"((A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409:0.00143,A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409:0.00070)NODE_00000161.00:0.01778,(A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416:0.00647,(A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409:0.00210,(A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423:0.00578,(A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409:0.00350,(A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412:0.00428,(A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409:0.00214,((A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409:0.00214,A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416:0.00357)NODE_00000090.98:0.00426,(A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428:0.00761,(A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416:0.00517,((A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409:0.00142,A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409:0.00214)NODE_00000050.93:0.00214,(A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409:0.00285,((A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409:0.00500,A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409:0.00071)NODE_00000010.47:0.00055,(A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409:0.00357,A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409:0.00142)NODE_00000020.45:0.00055)NODE_0000000:0.00071)NODE_00000030.81:0.00286)NODE_00000040.96:0.00486)NODE_00000060.98:0.00184)NODE_00000070.79:0.00199)NODE_00000080.81:0.00074)NODE_00000100.81:0.00071)NODE_00000110.80:0.00055)NODE_00000120.87:0.01374)NODE_00000131.00:0.01010)NODE_00000141.00:0.00066)NODE_00000150.73:0.00408)NODE_0000017:0.00100;";

  static INPUT_DATES: LazyLock<DatesMap> = LazyLock::new(|| {
    read_dates(
      "../../data/flu/h3n2/20/metadata.tsv",
      &Some(o!("name")),
      &Some(o!("date")),
    )
    .unwrap()
  });

  // Python TreeTime v0 baseline results obtained from:
  //
  // ./dev/docker/python \
  //   treetime \
  //   --branch-length-mode input \
  //   --clock-rate 0.002826 \
  //   --dates data/flu/h3n2/20/metadata.tsv \
  //   --keep-polytomies \
  //   --keep-root \
  //   --outdir tmp/python_treetime_baseline \
  //   --sequence-len 1407 \
  //   --tree tmp/clock_test/rerooted.newick \
  //
  // Output: tmp/python_treetime_baseline/dates.tsv
  static EXPECTED_TIMES_BRANCH_LEN_MODE_INPUT: LazyLock<BTreeMap<String, f64>> = LazyLock::new(|| {
    btreemap! {
      o!("NODE_0000017") => 1996.940903,
      o!("NODE_00000161.00") => 2002.844210,
      o!("NODE_00000150.73") => 1998.507899,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_00000141.00") => 1998.775604,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_00000131.00") => 2001.478893,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_00000120.87") => 2006.569842,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_00000110.80") => 2006.906449,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_00000100.81") => 2007.219376,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_00000080.81") => 2007.480821,
      o!("NODE_00000090.98") => 2008.633701,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_00000070.79") => 2008.485779,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_00000060.98") => 2009.180734,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_00000040.96") => 2010.631447,
      o!("NODE_00000050.93") => 2011.356438,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_00000030.81") => 2011.533989,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.877473,
      o!("NODE_00000010.47") => 2012.068368,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_00000020.45") => 2012.146432,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
    }
  });

  // Python TreeTime v0 baseline results obtained from:
  //
  // treetime \
  //   --tree tmp/clock_test/rerooted.newick \
  //   --keep-root \
  //   --clock-rate 0.002826 \
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
  static EXPECTED_TIMES_BRANCH_LEN_MODE_MARGINAL: LazyLock<BTreeMap<String, f64>> = LazyLock::new(|| {
    btreemap! {
      o!("NODE_0000017") => 1996.580809,
      o!("NODE_00000161.00") => 2002.837837,
      o!("A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409") => 2003.840520,
      o!("A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409") => 2003.002738,
      o!("NODE_00000150.73") => 1998.383572,
      o!("A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416") => 2000.681725,
      o!("NODE_00000141.00") => 1998.833856,
      o!("A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409") => 2000.134155,
      o!("NODE_00000131.00") => 2001.497951,
      o!("A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423") => 2003.002738,
      o!("NODE_00000120.87") => 2006.423551,
      o!("A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409") => 2007.487337,
      o!("NODE_00000110.80") => 2006.704928,
      o!("A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412") => 2008.865161,
      o!("NODE_00000100.81") => 2007.109002,
      o!("A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409") => 2008.150582,
      o!("NODE_00000080.81") => 2007.447819,
      o!("NODE_00000090.98") => 2008.645795,
      o!("A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409") => 2009.481862,
      o!("A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416") => 2009.522929,
      o!("NODE_00000070.79") => 2008.431561,
      o!("A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428") => 2011.980151,
      o!("NODE_00000060.98") => 2009.155837,
      o!("A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416") => 2011.651608,
      o!("NODE_00000040.96") => 2010.533103,
      o!("NODE_00000050.93") => 2011.358422,
      o!("A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409") => 2011.955510,
      o!("A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409") => 2011.985626,
      o!("NODE_00000030.81") => 2011.464074,
      o!("A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409") => 2012.257358,
      o!("NODE_0000000") => 2011.894667,
      o!("NODE_00000010.47") => 2012.088124,
      o!("A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409") => 2013.112252,
      o!("A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409") => 2012.856947,
      o!("NODE_00000020.45") => 2012.133696,
      o!("A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409") => 2013.405202,
      o!("A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409") => 2012.837782,
    }
  });

  static ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  static ALN: LazyLock<Vec<FastaRecord>> = LazyLock::new(|| {
    let alphabet = ALPHABET.clone();
    read_many_fasta(&["../../data/flu/h3n2/20/aln.fasta.xz"], &alphabet).unwrap()
  });

  #[test]
  fn test_timetree_flu_h3n2_poisson() -> Result<(), Report> {
    const CLOCK_RATE_MU: f64 = 0.0028;
    const SEQUENCE_LENGTH_L: usize = 1400;

    fn dump_graph(graph: &GraphAncestral, filename: &str) -> Result<(), Report> {
      let OUTPUT_DIR = Path::new("../../tmp/test_timetree_flu_h3n2_poisson");
      json_write_file(OUTPUT_DIR.join(filename), graph, JsonPretty(true))?;
      Ok(())
    }

    let graph: GraphAncestral = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;

    load_date_constraints(&INPUT_DATES, &graph)?;

    create_poisson_branch_distributions(&graph, CLOCK_RATE_MU, SEQUENCE_LENGTH_L, BRANCH_GRID_SIZE)?;
    dump_graph(&graph, "001_after_create_poisson_branch_distributions.json")?;

    propagate_distributions_backward(&graph)?;
    dump_graph(&graph, "002_after_propagate_distributions_backward.json")?;

    propagate_distributions_forward(&graph)?;
    dump_graph(&graph, "003_after_propagate_distributions_forward.json")?;

    let actual = round_times(&extract_node_times(&graph));
    let expected = round_times(&EXPECTED_TIMES_BRANCH_LEN_MODE_INPUT);

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
    fn dump_graph(graph: &GraphAncestral, filename: &str) -> Result<(), Report> {
      let OUTPUT_DIR = Path::new("../../tmp/test_timetree_flu_h3n2_marginal_sparse");
      json_write_file(OUTPUT_DIR.join(filename), graph, JsonPretty(true))?;
      Ok(())
    }

    let mut graph: GraphAncestral = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;

    load_date_constraints(&INPUT_DATES, &graph)?;

    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: SEQUENCE_LENGTH,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &ALN)?;
    dump_graph(&graph, "001_after_compress_sequences.json")?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll>>> = vec![sparse_partition];

    // Run marginal reconstruction
    run_marginal(&graph, &partitions, Some(&ALN))?;
    dump_graph(&graph, "002_after_run_marginal.json")?;

    initialize_node_divergences(&graph);
    dump_graph(&graph, "002_after_calculate_divs.json")?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockOptions::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
    )?;
    dump_graph(&graph, "003_after_clock_model.json")?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    dump_graph(&graph, "004_after_initialize_node_times.json")?;

    run_timetree(&mut graph, &partitions, true, &clock_model)?;
    dump_graph(&graph, "005_after_run_timetree.json")?;

    let actual = round_times(&extract_node_times(&graph));
    let expected = round_times(&EXPECTED_TIMES_BRANCH_LEN_MODE_MARGINAL);

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

    let mut graph: GraphAncestral = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;

    load_date_constraints(&INPUT_DATES, &graph)?;

    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: SEQUENCE_LENGTH_L,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll>>> = vec![dense_partition];

    run_marginal(&graph, &partitions, Some(&ALN))?;
    dump_graph(&graph, "001_after_run_marginal.json")?;

    initialize_node_divergences(&graph);
    dump_graph(&graph, "002_after_calculate_divs.json")?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockOptions::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
    )?;
    dump_graph(&graph, "003_after_clock_model.json")?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    dump_graph(&graph, "004_after_initialize_node_times.json")?;

    run_timetree(&mut graph, &partitions, true, &clock_model)?;
    dump_graph(&graph, "005_after_run_timetree.json")?;

    let actual = round_times(&extract_node_times(&graph));
    let expected = round_times(&EXPECTED_TIMES_BRANCH_LEN_MODE_MARGINAL);

    {
      assert_eq!(&expected.into_iter().collect_vec(), &actual.into_iter().collect_vec());
    }

    Ok(())
  }

  fn extract_node_times(graph: &GraphAncestral) -> BTreeMap<String, f64> {
    graph
      .get_nodes()
      .into_iter()
      .filter_map(|node_ref| {
        let node = node_ref.read_arc();
        let payload = node.payload().read_arc();
        let name = payload.name.as_ref()?.clone();
        let time = payload.time?;
        Some((name, time))
      })
      .collect()
  }

  fn round_to(x: f64, n: u32) -> f64 {
    let p = 10_f64.powi(n as i32);
    (x * p).round() / p
  }

  fn round_times(times: &BTreeMap<String, f64>) -> BTreeMap<String, f64> {
    times
      .iter()
      .map(|(name, time)| (name.clone(), round_to(*time, 4)))
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
