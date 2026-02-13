use clap::Parser;
use ctor::ctor;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::{Arc, LazyLock};
use treetime::alphabet::alphabet::Alphabet;
use treetime::commands::ancestral::fitch::compress_sequences;
use treetime::commands::ancestral::marginal::initialize_marginal;
use treetime::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use treetime::commands::clock::date_constraints::load_date_constraints;
use treetime::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use treetime::commands::timetree::inference::backward_pass::propagate_distributions_backward;
use treetime::commands::timetree::inference::forward_pass::propagate_distributions_forward;
use treetime::commands::timetree::inference::runner::{BRANCH_GRID_SIZE, run_timetree};
use treetime::commands::timetree::partition_ops::PartitionTimetreeAll;
use treetime::commands::timetree::utils::{
  initialize_clock_totals_from_time_distributions, initialize_node_divergences,
};
use treetime::distribution::distribution::Distribution;
use treetime::distribution::distribution_function::DistributionFunction;
use treetime::gtr::get_gtr::{JC69Params, jc69};
use treetime::io::dates_csv::{DatesMap, read_dates};
use treetime::io::fasta::{FastaRecord, read_many_fasta};
use treetime::io::nwk::nwk_read_str;
use treetime::o;
use treetime::representation::edge_timetree::EdgeTimetree;
use treetime::representation::node_timetree::NodeTimetree;
use treetime::representation::partition_marginal_dense::PartitionMarginalDense;
use treetime::representation::partition_marginal_sparse::PartitionMarginalSparse;
use treetime::representation::partition_timetree::GraphTimetree;
use treetime_graph::edge::HasBranchLength;
use treetime_graph::node::Named;
use treetime_io::json::{JsonPretty, json_write_file};
use treetime_utils::global_init::global_init;
use treetime_utils::string::truncate_right_with_ellipsis;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Clone, Serialize, Deserialize)]
#[command(
  name = "timetree-validation",
  about = "Validates TreeTime timetree inference against Python v0 baseline results",
  version
)]
pub struct Args {
  /// Tests to run (comma-separated: poisson, marginal_sparse, marginal_dense)
  #[arg(
    long,
    value_delimiter = ',',
    default_value = "poisson,marginal_sparse,marginal_dense"
  )]
  tests: Vec<String>,

  /// Output directory base path
  #[arg(long, default_value = "tmp")]
  output_dir: String,

  /// Absolute difference threshold for filtering comparison table rows
  #[arg(long, default_value_t = 0.001)]
  abs_diff_threshold: f64,

  /// Branch grid size for Poisson distribution discretization
  #[arg(long, default_value_t = BRANCH_GRID_SIZE)]
  branch_grid_size: usize,

  /// Optional statistics output file path
  #[arg(long)]
  stats_output: Option<String>,

  /// Enable detailed progress output
  #[arg(long, short = 'v')]
  verbose: bool,
}

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

static INPUT_DATES: LazyLock<DatesMap> =
  LazyLock::new(|| read_dates("data/flu/h3n2/20/metadata.tsv", &Some(o!("name")), &Some(o!("date"))).unwrap());

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
  read_many_fasta(&["data/flu/h3n2/20/aln.fasta.xz"], &alphabet).unwrap()
});

#[derive(Clone, Serialize, Deserialize)]
struct TestResult {
  name: String,
  expected: BTreeMap<String, f64>,
  actual: BTreeMap<String, f64>,
}

#[derive(Clone, Serialize, Deserialize)]
struct DiffStats {
  mean_abs_diff: f64,
  median_abs_diff: f64,
  max_abs_diff: f64,
  min_abs_diff: f64,
  std_dev: f64,
  rmse: f64,
}

fn main() -> Result<(), Report> {
  let args = Args::parse();
  let mut results = Vec::new();

  println!("# TreeTime Validation Results\n");

  if args.verbose {
    println!("Running tests: {}\n", args.tests.join(", "));
  }

  for test_name in &args.tests {
    match test_name.to_lowercase().as_str() {
      "poisson" => {
        let result = run_poisson_test(&args)?;
        print_comparison_table(&result, args.abs_diff_threshold);
        results.push(result);
      },
      "marginal_sparse" => {
        let result = run_marginal_sparse_test(&args)?;
        print_comparison_table(&result, args.abs_diff_threshold);
        results.push(result);
      },
      "marginal_dense" => {
        let result = run_marginal_dense_test(&args)?;
        print_comparison_table(&result, args.abs_diff_threshold);
        results.push(result);
      },
      _ => {
        eprintln!("Warning: Unknown test '{test_name}', skipping");
      },
    }
  }

  if results.is_empty() {
    eprintln!("Error: No valid tests were run");
    return Ok(());
  }

  println!("\n## Aggregate Statistics Summary\n");
  print_aggregate_statistics(&results);

  if let Some(stats_output_path) = &args.stats_output {
    write_statistics_to_file(stats_output_path, &results)?;
    if args.verbose {
      println!("Statistics written to: {stats_output_path}");
    }
  }

  Ok(())
}

fn dump_graph(graph: &GraphTimetree, output_dir: &str, filename: &str) -> Result<(), Report> {
  let output_path = Path::new(output_dir);
  json_write_file(output_path.join(filename), graph, JsonPretty(true))?;
  Ok(())
}

fn run_poisson_test(args: &Args) -> Result<TestResult, Report> {
  const CLOCK_RATE_MU: f64 = 0.0028;
  const SEQUENCE_LENGTH_L: usize = 1400;

  let output_dir_str = format!("{}/test_timetree_flu_h3n2_poisson", args.output_dir);

  println!("### Test 1: Poisson (branch length input mode)");

  let graph: GraphTimetree = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;
  load_date_constraints(&INPUT_DATES, &graph)?;

  create_poisson_branch_distributions(&graph, CLOCK_RATE_MU, SEQUENCE_LENGTH_L, args.branch_grid_size)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "001_after_create_poisson_branch_distributions.json",
  )?;

  propagate_distributions_backward(&graph, None)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "002_after_propagate_distributions_backward.json",
  )?;

  propagate_distributions_forward(&graph)?;
  dump_graph(
    &graph,
    &output_dir_str,
    "003_after_propagate_distributions_forward.json",
  )?;

  let actual = extract_node_times(&graph);
  let expected = EXPECTED_TIMES_BRANCH_LEN_MODE_INPUT.clone();

  println!("✓ Completed\n");

  Ok(TestResult {
    name: o!("Poisson (input mode)"),
    expected,
    actual,
  })
}

fn run_marginal_sparse_test(args: &Args) -> Result<TestResult, Report> {
  let output_dir_str = format!("{}/test_timetree_flu_h3n2_marginal_sparse", args.output_dir);

  println!("### Test 2: Marginal Sparse");

  let mut graph: GraphTimetree = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;
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
  dump_graph(&graph, &output_dir_str, "001_after_compress_sequences.json")?;

  let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![sparse_partition];

  initialize_marginal(&graph, &partitions, &ALN)?;
  dump_graph(&graph, &output_dir_str, "002_after_run_marginal.json")?;

  initialize_node_divergences(&graph);
  dump_graph(&graph, &output_dir_str, "002_after_calculate_divs.json")?;

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(CLOCK_RATE),
    true,
    &BranchPointOptimizationParams::default(),
  )?;
  dump_graph(&graph, &output_dir_str, "003_after_clock_model.json")?;

  initialize_clock_totals_from_time_distributions(&graph)?;
  dump_graph(&graph, &output_dir_str, "004_after_initialize_node_times.json")?;

  run_timetree(&mut graph, &partitions, &clock_model, None)?;
  dump_graph(&graph, &output_dir_str, "005_after_run_timetree.json")?;

  let actual = extract_node_times(&graph);
  let expected = EXPECTED_TIMES_BRANCH_LEN_MODE_MARGINAL.clone();

  println!("✓ Completed\n");

  Ok(TestResult {
    name: o!("Marginal Sparse"),
    expected,
    actual,
  })
}

fn run_marginal_dense_test(args: &Args) -> Result<TestResult, Report> {
  let output_dir_str = format!("{}/test_timetree_flu_h3n2_marginal_dense", args.output_dir);

  println!("### Test 3: Marginal Dense");

  let mut graph: GraphTimetree = nwk_read_str(INPUT_REROOTED_TREE_NWK)?;
  load_date_constraints(&INPUT_DATES, &graph)?;

  let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
    index: 0,
    gtr: jc69(JC69Params::default())?,
    alphabet: ALPHABET.clone(),
    length: SEQUENCE_LENGTH,
    nodes: btreemap! {},
    edges: btreemap! {},
  }));

  let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![dense_partition];

  initialize_marginal(&graph, &partitions, &ALN)?;
  dump_graph(&graph, &output_dir_str, "001_after_run_marginal.json")?;

  initialize_node_divergences(&graph);
  dump_graph(&graph, &output_dir_str, "002_after_calculate_divs.json")?;

  let clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &ClockParams::default(),
    Some(CLOCK_RATE),
    true,
    &BranchPointOptimizationParams::default(),
  )?;
  dump_graph(&graph, &output_dir_str, "003_after_clock_model.json")?;

  initialize_clock_totals_from_time_distributions(&graph)?;
  dump_graph(&graph, &output_dir_str, "004_after_initialize_node_times.json")?;

  run_timetree(&mut graph, &partitions, &clock_model, None)?;
  dump_graph(&graph, &output_dir_str, "005_after_run_timetree.json")?;

  let actual = extract_node_times(&graph);
  let expected = EXPECTED_TIMES_BRANCH_LEN_MODE_MARGINAL.clone();

  println!("✓ Completed\n");

  Ok(TestResult {
    name: o!("Marginal Dense"),
    expected,
    actual,
  })
}

fn extract_node_times(graph: &GraphTimetree) -> BTreeMap<String, f64> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let node = node_ref.read_arc();
      let payload = node.payload().read_arc();
      let name = payload.name()?.as_ref().to_owned();
      let time = payload.time?;
      Some((name, time))
    })
    .collect()
}

fn create_poisson_branch_distributions(
  graph: &GraphTimetree,
  mu: f64,
  seq_len: usize,
  n_points: usize,
) -> Result<(), Report> {
  let seq_len_f64 = seq_len as f64;

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();

    if let Some(branch_length) = edge.branch_length() {
      let expected_time = branch_length / mu;
      let max_time = 3.0 * expected_time.max(1.0);
      let dx = max_time / (n_points - 1) as f64;

      let y = Array1::from_shape_fn(n_points, |i| {
        let dt = i as f64 * dx;
        if dt < 1e-10 {
          0.0
        } else {
          let log_p = -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln();
          log_p.exp()
        }
      });

      let y_max = y.iter().copied().map(OrderedFloat).max().map_or(1.0, |x| x.0);
      let y_normalized = y.mapv(|v| v / y_max);

      let distribution_fn = DistributionFunction::from_start_dx_values(0.0, dx, y_normalized)?;
      let distribution = Distribution::Function(distribution_fn);
      edge.branch_length_distribution = Some(Arc::new(distribution));
    }
  }

  Ok(())
}

fn calculate_diff_stats(expected: &BTreeMap<String, f64>, actual: &BTreeMap<String, f64>) -> DiffStats {
  let diffs: Vec<f64> = expected
    .keys()
    .filter_map(|key| {
      let exp = expected.get(key)?;
      let act = actual.get(key)?;
      Some((exp - act).abs())
    })
    .collect();

  let n = diffs.len() as f64;
  let mean_abs_diff = diffs.iter().sum::<f64>() / n;

  let mut sorted_diffs = diffs.clone();
  sorted_diffs.sort_by_key(|&x| OrderedFloat(x));
  let median_abs_diff = if sorted_diffs.is_empty() {
    0.0
  } else {
    sorted_diffs[sorted_diffs.len() / 2]
  };

  let max_abs_diff = diffs.iter().copied().map(OrderedFloat).max().map_or(0.0, |x| x.0);
  let min_abs_diff = diffs.iter().copied().map(OrderedFloat).min().map_or(0.0, |x| x.0);

  let variance = diffs.iter().map(|&d| (d - mean_abs_diff).powi(2)).sum::<f64>() / n;
  let std_dev = variance.sqrt();

  let squared_diffs: f64 = expected
    .keys()
    .filter_map(|key| {
      let exp = expected.get(key)?;
      let act = actual.get(key)?;
      Some((exp - act).powi(2))
    })
    .sum();
  let rmse = (squared_diffs / n).sqrt();

  DiffStats {
    mean_abs_diff,
    median_abs_diff,
    max_abs_diff,
    min_abs_diff,
    std_dev,
    rmse,
  }
}

fn print_comparison_table(result: &TestResult, abs_diff_threshold: f64) {
  println!("### {}\n", result.name);

  let stats = calculate_diff_stats(&result.expected, &result.actual);

  let mut rows: Vec<_> = result
    .expected
    .keys()
    .map(|key| {
      let expected = result.expected.get(key).copied().unwrap_or(0.0);
      let actual = result.actual.get(key).copied().unwrap_or(0.0);
      let diff = (expected - actual).abs();
      (key.clone(), expected, actual, diff)
    })
    .collect();

  rows.sort_by_key(|(_, _, _, diff)| OrderedFloat(-diff));

  let filtered_rows: Vec<_> = rows
    .iter()
    .filter(|(_, _, _, diff)| *diff >= abs_diff_threshold)
    .collect();
  let hidden_count = rows.len() - filtered_rows.len();

  println!(
    "| {:<50} | {:>12} | {:>12} | {:>12} |",
    "Node", "Expected", "Actual", "Abs Diff"
  );
  println!("| {:-<50} | {:-<12} | {:-<12} | {:-<12} |", "", "", "", "");

  for (node, expected, actual, diff) in &filtered_rows {
    let short_name = truncate_right_with_ellipsis(node, 50);
    println!("| {short_name:<50} | {expected:>12.6} | {actual:>12.6} | {diff:>12.6} |");
  }

  if hidden_count > 0 {
    println!();
    println!("({hidden_count} nodes with abs diff < {abs_diff_threshold:.6} are hidden)");
  }

  println!();
  println!("#### Statistics");
  println!();
  println!("| {:<30} | {:>12} |", "Metric", "Value");
  println!("| {:-<30} | {:-<12} |", "", "");
  println!("| {:<30} | {:>12.6} |", "Mean Absolute Difference", stats.mean_abs_diff);
  println!(
    "| {:<30} | {:>12.6} |",
    "Median Absolute Difference", stats.median_abs_diff
  );
  println!("| {:<30} | {:>12.6} |", "Max Absolute Difference", stats.max_abs_diff);
  println!("| {:<30} | {:>12.6} |", "Min Absolute Difference", stats.min_abs_diff);
  println!("| {:<30} | {:>12.6} |", "Standard Deviation", stats.std_dev);
  println!("| {:<30} | {:>12.6} |", "RMSE", stats.rmse);
  println!();
}

fn print_aggregate_statistics(results: &[TestResult]) {
  println!(
    "| {:<25} | {:>15} | {:>17} | {:>14} | {:>10} | {:>11} |",
    "Test", "Mean Abs Diff", "Median Abs Diff", "Max Abs Diff", "RMSE", "Std Dev"
  );
  println!(
    "| {:-<25} | {:-<15} | {:-<17} | {:-<14} | {:-<10} | {:-<11} |",
    "", "", "", "", "", ""
  );

  for result in results {
    let stats = calculate_diff_stats(&result.expected, &result.actual);
    println!(
      "| {:<25} | {:>15.6} | {:>17.6} | {:>14.6} | {:>10.6} | {:>11.6} |",
      result.name, stats.mean_abs_diff, stats.median_abs_diff, stats.max_abs_diff, stats.rmse, stats.std_dev
    );
  }

  let overall_mean = results
    .iter()
    .map(|r| calculate_diff_stats(&r.expected, &r.actual).mean_abs_diff)
    .sum::<f64>()
    / results.len() as f64;

  let overall_rmse = results
    .iter()
    .map(|r| calculate_diff_stats(&r.expected, &r.actual).rmse)
    .sum::<f64>()
    / results.len() as f64;

  println!(
    "| {:-<25} | {:-<15} | {:-<17} | {:-<14} | {:-<10} | {:-<11} |",
    "", "", "", "", "", ""
  );
  println!(
    "| {:<25} | {:>15.6} | {:>17} | {:>14} | {:>10.6} | {:>11} |",
    "Overall", overall_mean, "-", "-", overall_rmse, "-"
  );
  println!();
}

fn write_statistics_to_file(path: &str, results: &[TestResult]) -> Result<(), Report> {
  use std::fs::File;
  use std::io::Write;

  let mut file = File::create(path)?;

  writeln!(
    file,
    "Test,Mean Abs Diff,Median Abs Diff,Max Abs Diff,Min Abs Diff,Std Dev,RMSE"
  )?;

  for result in results {
    let stats = calculate_diff_stats(&result.expected, &result.actual);
    writeln!(
      file,
      "{},{},{},{},{},{},{}",
      result.name,
      stats.mean_abs_diff,
      stats.median_abs_diff,
      stats.max_abs_diff,
      stats.min_abs_diff,
      stats.std_dev,
      stats.rmse
    )?;
  }

  Ok(())
}
