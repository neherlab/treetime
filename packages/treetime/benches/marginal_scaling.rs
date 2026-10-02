#![allow(
  clippy::expect_used,
  clippy::unwrap_used,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ctor::ctor;
use rayon::ThreadPoolBuilder;
use std::hint::black_box;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::mask::create_mask;
use treetime::ancestral::params::{AncestralParams, MethodAncestral};
use treetime::ancestral::pipeline::run;
use treetime::ancestral::sample::SampleMode;
use treetime::cancel::NoopCancel;
use treetime::gtr::get_gtr::GtrModelName;
use treetime::progress::NoopProgress;
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, get_common_length, node_seq_inputs};
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::AlignmentRecord;
use treetime_utils::init::global::global_init;

const DATASET_SEQUENCES: u64 = 200;

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn benchmark_marginal_scaling(criterion: &mut Criterion) {
  let mut group = criterion.benchmark_group("marginal_reconstruction_threads");
  group.sample_size(10);
  group.throughput(Throughput::Elements(DATASET_SEQUENCES));

  let (input, mask) = setup();
  let alphabet = Alphabet::default();
  let params = AncestralParams {
    method: MethodAncestral::Marginal,
    model: GtrModelName::JC69,
    dense: Some(false),
    include_leaves: false,
    impute_missing_data: false,
    gtr_iterations: 0,
    site_specific_gtr: false,
    seed: Some(0),
    sample_from_profile: SampleMode::Argmax,
  };

  for threads in [1, 2, 4, 8] {
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    group.bench_with_input(BenchmarkId::new("sparse", threads), &threads, |bencher, _| {
      bencher.iter(|| {
        let output = pool
          .install(|| {
            run(
              &params,
              black_box(&input),
              alphabet.clone(),
              mask.clone(),
              &NoopCancel,
              &NoopProgress,
              &NoopProgress,
            )
          })
          .unwrap();
        black_box(output);
      });
    });
  }
  group.finish();
}

fn setup() -> (AncestralInput, Vec<bool>) {
  let alphabet = Alphabet::default();
  let project_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
  let nwk_parsed = nwk_read_file(project_root.join("data/flu/h3n2/200/tree.nwk")).unwrap();
  let names = nwk_parsed.names();
  let alignment: Vec<AlignmentRecord> =
    read_many_fasta_path(&[project_root.join("data/flu/h3n2/200/aln.fasta.xz")], &alphabet)
      .unwrap()
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
  let mask = create_mask(&alignment, get_common_length(&alignment).unwrap(), &alphabet);
  let input = AncestralInput {
    nodes: node_seq_inputs(&nwk_parsed.graph, &names, alignment),
    edges: nwk_parsed
      .branch_lengths
      .into_iter()
      .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
      .collect(),
    graph: nwk_parsed.graph,
  };
  (input, mask)
}

criterion_group!(benches, benchmark_marginal_scaling);
criterion_main!(benches);
