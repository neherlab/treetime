#![allow(
  clippy::expect_used,
  clippy::unwrap_used,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

use criterion::{BatchSize, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ctor::ctor;
use eyre::Report;
use rayon::ThreadPoolBuilder;
use std::hint::black_box;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::params::{AncestralParams, MethodAncestral};
use treetime::ancestral::pipeline::run;
use treetime::cancel::NoopCancel;
use treetime::gtr::get_gtr::GtrModelName;
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::NoopProgress;
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, pair_leaf_sequences};
use treetime::seq::sink::{SeqItem, SeqSink};
use treetime_io::fasta::fasta_read_file;
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

  let params = AncestralParams {
    method: MethodAncestral::Marginal,
    model: GtrModelName::JC69,
    dense: Some(false),
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: false,
    gtr_iterations: 0,
    site_specific_gtr: false,
    seed: 0,
    sample_from_profile: SampleMode::Argmax,
  };

  for threads in [1, 2, 4, 8] {
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    group.bench_with_input(BenchmarkId::new("sparse", threads), &threads, |bencher, _| {
      bencher.iter_batched(
        setup,
        |input| {
          let output = pool
            .install(|| {
              run(
                &params,
                black_box(input),
                Some(&mut DiscardSequences),
                &NoopCancel,
                &NoopProgress,
                &NoopProgress,
              )
            })
            .unwrap();
          black_box(output);
        },
        BatchSize::LargeInput,
      );
    });
  }
  group.finish();
}

struct DiscardSequences;

impl SeqSink for DiscardSequences {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    black_box(item.seq);
    Ok(())
  }
}

fn setup() -> AncestralInput {
  let alphabet = Alphabet::default();
  let project_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
  let nwk_parsed = nwk_read_file(project_root.join("data/flu/h3n2/200/tree.nwk")).unwrap();
  let names = nwk_parsed.names();
  let alignment: Vec<AlignmentRecord> = fasta_read_file(project_root.join("data/flu/h3n2/200/aln.fasta.xz"), &alphabet)
    .unwrap()
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
  let sequences = pair_leaf_sequences(&nwk_parsed.graph, &names, alignment).sequences;
  let mask = sequences.mask(sequences.common_length().unwrap(), &alphabet);
  AncestralInput {
    nodes: sequences.nodes,
    edges: nwk_parsed
      .branch_lengths
      .into_iter()
      .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
      .collect(),
    graph: nwk_parsed.graph,
    alphabet,
    mask,
  }
}

criterion_group!(benches, benchmark_marginal_scaling);
criterion_main!(benches);
