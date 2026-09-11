use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ctor::ctor;
use rayon::ThreadPoolBuilder;
use std::collections::BTreeMap;
use std::hint::black_box;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::fitch::create_fitch_partition;
use treetime::ancestral::marginal::{marginal_update, profile_branch_lengths};
use treetime::gtr::get_gtr::{JC69Params, jc69};
use treetime::payload::ancestral::GraphAncestral;
use treetime_graph::edge::GraphEdgeKey;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nwk::{NwkParse, nwk_read_file};
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
    let (graph, mut partitions, branch_lengths) = setup();
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    group.bench_with_input(BenchmarkId::new("sparse", threads), &threads, |bencher, _| {
      bencher.iter(|| {
        pool
          .install(|| {
            marginal_update(
              black_box(&graph),
              &profile_branch_lengths(black_box(&branch_lengths)),
              black_box(&mut partitions),
            )
          })
          .unwrap();
      });
    });
  }
  group.finish();
}

fn setup() -> (
  GraphAncestral,
  [treetime::partition::marginal::sparse::partition::PartitionMarginalSparse; 1],
  BTreeMap<GraphEdgeKey, Option<f64>>,
) {
  ThreadPoolBuilder::new()
    .num_threads(1)
    .build()
    .unwrap()
    .install(setup_inner)
}

fn setup_inner() -> (
  GraphAncestral,
  [treetime::partition::marginal::sparse::partition::PartitionMarginalSparse; 1],
  BTreeMap<GraphEdgeKey, Option<f64>>,
) {
  let alphabet = Alphabet::default();
  let project_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = nwk_read_file(project_root.join("data/flu/h3n2/200/tree.nwk")).unwrap();
  let alignment = read_many_fasta(&[project_root.join("data/flu/h3n2/200/aln.fasta.xz")], &alphabet).unwrap();
  let fitch = create_fitch_partition(&graph, 0, alphabet, &alignment, &names).unwrap();
  let gtr = jc69(JC69Params::default()).unwrap();
  let partition = fitch.into_marginal_sparse(gtr, &graph).unwrap();
  let mut partitions = [partition];
  marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions).unwrap();
  (graph, partitions, branch_lengths)
}

criterion_group!(benches, benchmark_marginal_scaling);
criterion_main!(benches);
