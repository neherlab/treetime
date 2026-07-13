use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use ctor::ctor;
use rayon::ThreadPoolBuilder;
use std::hint::black_box;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::fitch::create_fitch_partition;
use treetime::payload::ancestral::GraphAncestral;
use treetime_io::fasta::FastaRecord;
use treetime_io::nwk::nwk_read_str;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::init::global::global_init;

#[ctor]
fn init() {
  global_init();
}

fn benchmark_fitch_scaling(criterion: &mut Criterion) {
  let mut group = criterion.benchmark_group("fitch_informative_positions");
  group.sample_size(10);
  let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();

  for topology in [Topology::Balanced, Topology::Caterpillar, Topology::Star] {
    for length in [1_000, 10_000] {
      for (density, denominator) in [("invariant", None), ("one-percent", Some(100)), ("all", Some(1))] {
        let (graph, alignment) = setup(topology, 16, length, denominator);
        let benchmark = format!("{}-{density}", topology.name());
        group.throughput(Throughput::Elements((alignment.len() * length) as u64));
        group.bench_with_input(BenchmarkId::new(benchmark, length), &length, |bencher, _| {
          bencher.iter(|| {
            pool
              .install(|| create_fitch_partition(black_box(&graph), 0, Alphabet::default(), black_box(&alignment)))
              .unwrap();
          });
        });
      }
    }
  }
  group.finish();
}

#[derive(Clone, Copy)]
enum Topology {
  Balanced,
  Caterpillar,
  Star,
}

impl Topology {
  fn name(self) -> &'static str {
    match self {
      Self::Balanced => "balanced",
      Self::Caterpillar => "caterpillar",
      Self::Star => "star",
    }
  }
}

fn setup(
  topology: Topology,
  tips: usize,
  length: usize,
  informative_denominator: Option<usize>,
) -> (GraphAncestral, Vec<FastaRecord>) {
  let names = (0..tips).map(|index| format!("T{index}")).collect::<Vec<_>>();
  let newick = format!("{}root:0.0;", topology_newick(topology, &names));
  let graph = nwk_read_str(&newick).unwrap();
  let alignment = names
    .into_iter()
    .enumerate()
    .map(|(taxon, seq_name)| FastaRecord {
      seq_name,
      seq: (0..length)
        .map(|pos| {
          if informative_denominator.is_some_and(|denominator| pos % denominator == 0) && taxon % 2 == 1 {
            AsciiChar::from_byte_unchecked(b'C')
          } else {
            AsciiChar::from_byte_unchecked(b'A')
          }
        })
        .collect::<Seq>(),
      ..FastaRecord::default()
    })
    .collect();
  (graph, alignment)
}

fn topology_newick(topology: Topology, names: &[String]) -> String {
  match topology {
    Topology::Balanced => balanced_newick(names),
    Topology::Caterpillar => caterpillar_newick(names),
    Topology::Star => format!(
      "({})",
      names
        .iter()
        .map(|name| format!("{name}:0.1"))
        .collect::<Vec<_>>()
        .join(",")
    ),
  }
}

fn balanced_newick(names: &[String]) -> String {
  if names.len() == 1 {
    return names[0].clone();
  }
  let middle = names.len() / 2;
  format!(
    "({}:0.1,{}:0.1)",
    balanced_newick(&names[..middle]),
    balanced_newick(&names[middle..])
  )
}

fn caterpillar_newick(names: &[String]) -> String {
  match names {
    [name] => name.clone(),
    [first, rest @ ..] => format!("({first}:0.1,{}:0.1)", caterpillar_newick(rest)),
    [] => unreachable!(),
  }
}

criterion_group!(benches, benchmark_fitch_scaling);
criterion_main!(benches);
