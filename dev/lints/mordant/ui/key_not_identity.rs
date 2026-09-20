// Bad cases mirror real Scarlet compiler bugs: f64 to_bits as a dedup key,
// pointer bits as identity, and a span type as a map key. The test dylint.toml
// denies `Span` and enables both expression forms; a default (unconfigured)
// run stays silent on all of this.

use std::collections::{HashMap, HashSet};

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
struct Span {
    lo: u32,
    hi: u32,
}

fn span_key_is_flagged(spans: &[Span]) -> HashMap<Span, usize> {
    let mut m = HashMap::new();
    for (i, s) in spans.iter().enumerate() {
        m.insert(*s, i);
    }
    m
}

fn span_key_lookup_is_flagged(m: &HashMap<Span, usize>, s: Span) -> Option<usize> {
    m.get(&s).copied()
}

fn to_bits_key_is_flagged(constants: &[f64]) -> HashMap<u64, usize> {
    let mut pool = HashMap::new();
    for (i, c) in constants.iter().enumerate() {
        pool.insert(c.to_bits(), i);
    }
    pool
}

fn f32_to_bits_key_is_flagged(x: f32) -> HashSet<u32> {
    let mut seen = HashSet::new();
    seen.insert(x.to_bits());
    seen
}

fn pointer_cast_key_is_flagged(values: &[Box<str>]) -> HashMap<usize, ()> {
    let mut ids = HashMap::new();
    for v in values {
        ids.insert(v.as_ptr() as usize, ());
    }
    ids
}

struct Value(u64);

impl Value {
    fn to_raw(&self) -> u64 {
        self.0
    }
}

fn denied_method_key_is_flagged(vs: &[Value]) -> HashMap<u64, usize> {
    let mut m = HashMap::new();
    for (i, v) in vs.iter().enumerate() {
        m.insert(v.to_raw(), i);
    }
    m
}

fn method_as_value_is_fine(v: Value) -> HashMap<usize, u64> {
    let mut m = HashMap::new();
    m.insert(0, v.to_raw());
    m
}

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
struct FileId(u32);

fn composite_key_is_flagged(spans: &[Span]) -> HashMap<(Span, u32), usize> {
    let mut m = HashMap::new();
    for (i, s) in spans.iter().enumerate() {
        m.insert((*s, 0u32), i);
    }
    m
}

fn fixed_composite_is_fine(spans: &[Span]) -> HashMap<(FileId, Span), usize> {
    let mut m = HashMap::new();
    for (i, s) in spans.iter().enumerate() {
        m.insert((FileId(0), *s), i);
    }
    m
}

fn plain_u64_key_is_fine(lens: &[u64]) -> HashMap<u64, usize> {
    let mut m = HashMap::new();
    for (i, len) in lens.iter().enumerate() {
        m.insert(*len, i);
    }
    m
}

fn to_bits_as_value_is_fine(x: f64) -> HashMap<usize, u64> {
    let mut m = HashMap::new();
    m.insert(0, x.to_bits());
    m
}

fn span_in_value_position_is_fine(spans: &[Span]) -> HashMap<usize, Span> {
    let mut m = HashMap::new();
    for (i, s) in spans.iter().enumerate() {
        m.insert(i, *s);
    }
    m
}

fn main() {
    let spans = [Span { lo: 0, hi: 1 }];
    let m = span_key_is_flagged(&spans);
    let _ = span_key_lookup_is_flagged(&m, spans[0]);
    let _ = to_bits_key_is_flagged(&[1.0, 2.0]);
    let _ = f32_to_bits_key_is_flagged(1.0);
    let _ = pointer_cast_key_is_flagged(&[]);
    let _ = denied_method_key_is_flagged(&[Value(1)]);
    let _ = composite_key_is_flagged(&spans);
    let _ = fixed_composite_is_fine(&spans);
    let _ = method_as_value_is_fine(Value(2));
    let _ = plain_u64_key_is_fine(&[1, 2]);
    let _ = to_bits_as_value_is_fine(1.0);
    let _ = span_in_value_position_is_fine(&spans);
}
