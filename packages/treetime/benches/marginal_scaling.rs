use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ctor::ctor;
use rayon::ThreadPoolBuilder;
use std::collections::BTreeMap;
use std::hint::black_box;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::fitch::create_fitch_partition;
use treetime::ancestral::marginal::profile_branch_lengths;
use treetime::ancestral::pipeline::SparseReconstruction;
use treetime::gtr::get_gtr::{JC69Params, jc69};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::nwk::nwk_fasta_node_inputs;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::init::global::global_init;

#[ctor]
fn init() {
  global_init();
}

fn benchmark_marginal_scaling(criterion: &mut Criterion) {
  let mut group = criterion.benchmark_group("marginal_update_threads");
  group.sample_size(10);
  group.throughput(Throughput::Elements(200));

  for threads in [1, 2, 4, 8] {
    let (graph, recon, branch_lengths) = setup();
    // The update consumes the reconstruction and returns the next one, so the bench carries it from
    // iteration to iteration in a slot rather than updating one value in place.
    let mut slot = Some(recon);
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    group.bench_with_input(BenchmarkId::new("sparse", threads), &threads, |bencher, _| {
      bencher.iter(|| {
        let recon = slot
          .take()
          .expect("reconstruction is present at the start of an iteration");
        let (recon, _) = pool
          .install(|| recon.marginal_update(black_box(&graph), &profile_branch_lengths(black_box(&branch_lengths))))
          .unwrap();
        slot = Some(black_box(recon));
      });
    });
  }
  group.finish();
}

fn setup() -> (Graph, SparseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>) {
  ThreadPoolBuilder::new()
    .num_threads(1)
    .build()
    .unwrap()
    .install(setup_inner)
}

fn setup_inner() -> (Graph, SparseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>) {
  let alphabet = Alphabet::default();
  let project_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
  let nwk_parsed = nwk_read_file(project_root.join("data/flu/h3n2/200/tree.nwk")).unwrap();
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let alignment = read_many_fasta_path(&[project_root.join("data/flu/h3n2/200/aln.fasta.xz")], &alphabet).unwrap();
  let fitch = create_fitch_partition(&graph, 0, alphabet, &nwk_fasta_node_inputs(&graph, &names, alignment)).unwrap();
  let gtr = jc69(JC69Params::default()).unwrap();
  let (partition, node_states) = fitch.into_marginal_sparse(&graph).unwrap();
  let recon = SparseReconstruction::seeded(partition, gtr, node_states);
  let (recon, _) = recon
    .marginal_update(&graph, &profile_branch_lengths(&branch_lengths))
    .unwrap();
  (graph, recon, branch_lengths)
}

criterion_group!(benches, benchmark_marginal_scaling);
criterion_main!(benches);
