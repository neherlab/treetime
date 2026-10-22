#![expect(
  clippy::unwrap_used,
  reason = "benchmark setup: a failed step aborts the benchmark, as it fails a test"
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
use treetime::gtr::get_gtr::{GtrModel, JC69Params};
use treetime::partition::marginal::sample::SampleMode;
use treetime::progress::NoopProgress;
use treetime::seq::alignment::{AncestralInput, EdgeSeqInput, SequenceSource, pair_leaf_sequences};
use treetime::seq::sink::{SeqItem, SeqSink, SeqStream};
use treetime_graph::node_names::NodeNames;
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
    model: GtrModel::JC69(JC69Params::default()),
    dense: Some(false),
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: false,
    seed: 0,
    sample_from_profile: SampleMode::Argmax,
  };

  for threads in [1, 2, 4, 8] {
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    group.bench_with_input(BenchmarkId::new("sparse", threads), &threads, |bencher, _| {
      bencher.iter_batched(
        setup,
        |(input, names)| {
          let output = pool
            .install(|| {
              run(
                &params,
                black_box(input),
                Some(SeqStream::new(&mut DiscardSequences, &names)),
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

fn setup() -> (AncestralInput, NodeNames) {
  let alphabet = Alphabet::default();
  let project_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
  let nwk_parsed = nwk_read_file(project_root.join("data/flu/h3n2/200/tree.nwk")).unwrap();
  let names = nwk_parsed.names;
  let alignment: Vec<AlignmentRecord> = fasta_read_file(project_root.join("data/flu/h3n2/200/aln.fasta.xz"), &alphabet)
    .unwrap()
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
  let sequences = pair_leaf_sequences(&nwk_parsed.graph, &names, alignment).sequences;
  let length = sequences.common_length().unwrap();
  let mask = sequences.mask(length, &alphabet);
  let input = AncestralInput {
    edges: nwk_parsed
      .branch_lengths
      .into_iter()
      .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
      .collect(),
    graph: nwk_parsed.graph,
    sources: vec![SequenceSource {
      path: project_root.join("data/flu/h3n2/200/aln.fasta.xz"),
      alphabet,
      length,
      nodes: sequences.nodes,
      mask,
    }],
  };
  (input, names)
}

criterion_group!(benches, benchmark_marginal_scaling);
criterion_main!(benches);
