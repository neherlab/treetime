use criterion::{black_box, criterion_group, criterion_main, Criterion};
use treetime::utils::random::random_sequence;

fn vec_to_string_collect(v: Vec<char>) -> String {
  v.into_iter().collect()
}

fn vec_to_string_extend(v: Vec<char>) -> String {
  let mut s = String::new();
  s.extend(v);
  s
}

fn vec_to_string_from_iter(v: Vec<char>) -> String {
  String::from_iter(v)
}

fn vec_to_string_from_utf8(v: Vec<char>) -> String {
  let bytes: Vec<u8> = v.into_iter().map(|c| c as u8).collect();
  String::from_utf8(bytes).unwrap()
}

fn vec_to_string_push(v: Vec<char>) -> String {
  let mut s = String::with_capacity(v.len());
  for c in v {
    s.push(c);
  }
  s
}

fn benchmark_main(c: &mut Criterion) {
  let seq = random_sequence(1_000_000, Some(42));

  let mut g = c.benchmark_group("vec_to_string");

  g.bench_function("vec_to_string_collect", |b| {
    b.iter(|| vec_to_string_collect(black_box(seq.clone())));
  });

  g.bench_function("vec_to_string_extend", |b| {
    b.iter(|| vec_to_string_extend(black_box(seq.clone())));
  });

  g.bench_function("vec_to_string_from_iter", |b| {
    b.iter(|| vec_to_string_from_iter(black_box(seq.clone())));
  });

  g.bench_function("vec_to_string_from_utf8", |b| {
    b.iter(|| vec_to_string_from_utf8(black_box(seq.clone())));
  });

  g.bench_function("vec_to_string_push", |b| {
    b.iter(|| vec_to_string_push(black_box(seq.clone())));
  });

  g.finish();
}

criterion_group!(vec_char_to_string, benchmark_main);
criterion_main!(vec_char_to_string);
