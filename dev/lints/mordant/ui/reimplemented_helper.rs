// A function whose signature and body repeat another function's is one
// helper written twice; anything that differs in type or shape is not.

pub fn clamp_add(base: u32, delta: u32, limit: u32) -> u32 {
    let sum = base.saturating_add(delta);
    if sum > limit { limit } else { sum }
}

// Flagged: `clamp_add` again with the parameters renamed.
pub fn bounded_sum(a: u32, b: u32, max: u32) -> u32 {
    let total = a.saturating_add(b);
    if total > max { max } else { total }
}

// Fine: same shape, but the parameters are `u64`, so the signature differs.
pub fn clamp_add_wide(base: u64, delta: u64, limit: u64) -> u64 {
    let sum = base.saturating_add(delta);
    if sum > limit { limit } else { sum }
}

// Fine: one operator differs.
pub fn clamp_add_inclusive(base: u32, delta: u32, limit: u32) -> u32 {
    let sum = base.saturating_add(delta);
    if sum >= limit { limit } else { sum }
}

pub struct Row {
    cells: Vec<u32>,
    pad: u32,
}

pub struct Column {
    cells: Vec<u32>,
    pad: u32,
}

impl Row {
    pub fn extent(&self) -> u32 {
        self.cells.iter().sum::<u32>() + self.pad * 2 + 1
    }

    // Flagged: the same computation as `extent` on the same type.
    pub fn footprint(&self) -> u32 {
        self.cells.iter().sum::<u32>() + self.pad * 2 + 1
    }

    // Fine: too small to compare (under `reimplemented-helper-min-nodes`).
    pub fn pad(&self) -> u32 {
        self.pad
    }

    // Fine: an accessor as small as `pad`.
    pub fn margin(&self) -> u32 {
        self.pad
    }
}

impl Column {
    // Fine: spelled like `Row::extent`, but `self` is a `Column`, so the
    // signatures differ and the field reads are of another type.
    pub fn extent(&self) -> u32 {
        self.cells.iter().sum::<u32>() + self.pad * 2 + 1
    }
}

pub trait Resize {
    fn grow(&self, from: u32, to: u32) -> u32;
    fn shrink(&self, from: u32, to: u32) -> u32;
}

impl Row {
    fn remap(&self, from: u32, to: u32) -> u32 {
        if to > from {
            self.pad + (to - from)
        } else {
            self.pad.saturating_sub(from - to)
        }
    }
}

// Fine: the trait obliges the impl to define both, and each already forwards
// to the one shared helper, so there is no copy to delete.
impl Resize for Row {
    fn grow(&self, from: u32, to: u32) -> u32 {
        self.remap(from.min(to), to.max(from)) + self.cells.len() as u32
    }
    fn shrink(&self, from: u32, to: u32) -> u32 {
        self.remap(from.min(to), to.max(from)) + self.cells.len() as u32
    }
}

pub fn low_third((lo, _): (u32, u32), bias: u32) -> u32 {
    let v = lo.saturating_mul(3).saturating_add(bias);
    if v > 255 { 255 } else { v }
}

// Fine: the same body over the other half of the pair. The parameter
// patterns differ, so `hi` is not `lo` renamed.
pub fn high_third((_, hi): (u32, u32), bias: u32) -> u32 {
    let v = hi.saturating_mul(3).saturating_add(bias);
    if v > 255 { 255 } else { v }
}

// Flagged: destructures the same half as `low_third`.
pub fn first_third((first, _): (u32, u32), k: u32) -> u32 {
    let n = first.saturating_mul(3).saturating_add(k);
    if n > 255 { 255 } else { n }
}

pub struct Extent {
    pub start: u32,
    pub end: u32,
}

pub fn lead(Extent { start, .. }: Extent, pad: u32) -> u32 {
    start.saturating_sub(pad).saturating_mul(2).min(4096) / 8 + pad * 2
}

// Fine: reads the other field of `Extent` through a pattern of the same
// shape, so `end` is not `start` renamed.
pub fn trail(Extent { end, .. }: Extent, pad: u32) -> u32 {
    end.saturating_sub(pad).saturating_mul(2).min(4096) / 8 + pad * 2
}

pub trait Near {
    fn step(&self) -> u32;
}

pub trait Far {
    fn step(&self) -> u32;
}

pub fn walk_near<T: Near>(t: &T, from: u32, budget: u32) -> u32 {
    let hop = t.step().max(1);
    from.saturating_add(hop.saturating_mul(budget)) / hop + 1
}

// Fine: `t.step()` is `Far::step` here and `Near::step` above; the bound is
// part of the signature, so the two are different functions spelled alike.
pub fn walk_far<T: Far>(t: &T, from: u32, budget: u32) -> u32 {
    let hop = t.step().max(1);
    from.saturating_add(hop.saturating_mul(budget)) / hop + 1
}

// Quiet: closures are never compared structurally, so two bodies built
// around one never pair up even when they are copies.
pub fn doubled_evens(xs: &[u32]) -> Vec<u32> {
    xs.iter().filter(|x| **x % 2 == 0).map(|x| x * 2).collect()
}

// Quiet: the copy of `doubled_evens`.
pub fn twice_the_evens(values: &[u32]) -> Vec<u32> {
    values
        .iter()
        .filter(|v| **v % 2 == 0)
        .map(|v| v * 2)
        .collect()
}

fn main() {}
