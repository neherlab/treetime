#![allow(clippy::unwrap_used, reason = "benchmark code: a failed read aborts the benchmark")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::fmt::Write;
use std::hint::black_box;
use util_newick::{NewickDialect, NewickReadOptions, newick_from_str};

const LEAF_COUNT: usize = 20_000;

fn bench_read_dialects(c: &mut Criterion) {
  let mut group = c.benchmark_group("read_annotated_tree");
  group.sample_size(10);
  let inputs = [
    ("beast_tree", NewickDialect::BEAST, annotated_ladder(LEAF_COUNT, "")),
    (
      "enewick_beast_network",
      NewickDialect::ENEWICK_BEAST,
      annotated_ladder(LEAF_COUNT, "#H1"),
    ),
  ];
  for (label, dialect, input) in inputs {
    let options = NewickReadOptions {
      dialect,
      ..NewickReadOptions::default()
    };
    group.throughput(Throughput::Bytes(u64::try_from(input.len()).unwrap()));
    group.bench_with_input(BenchmarkId::from_parameter(label), &input, |bencher, input| {
      bencher.iter(|| newick_from_str(black_box(input), &options).unwrap());
    });
  }
  group.finish();
}

fn annotated_ladder(leaf_count: usize, hybrid: &str) -> String {
  let mut text = "(".repeat(leaf_count);
  write!(text, "(L0[&rate=1.25]:0.1){hybrid}[&segments={{0,1}}]:0.1").unwrap();
  for i in 1..leaf_count {
    write!(
      text,
      ",L{i}[&rate=1.25,height=0.5,height_95%_HPD={{0.25,0.75}}]:0.1)[&posterior=0.99,rate=1.0]:[&length=1.5]0.01"
    )
    .unwrap();
  }
  write!(text, ",{hybrid}[&segments={{1}}]:0.2);").unwrap();
  text
}

criterion_group!(benches, bench_read_dialects);
criterion_main!(benches);
