// A private fn returning a tuple with same-typed members that every caller
// destructures under the same names: the names exist everywhere except in
// the type, which is the one place that would stop a transposition.

struct Node {
    depth: u32,
    width: u32,
}

// Flagged: both callers, in two functions, bind the members `(depth, width)`;
// both are `u32`, so nothing but the names keeps them apart.
fn measure(n: &Node) -> (u32, u32) {
    (n.depth + 1, n.width * 2)
}

fn taller(a: &Node, b: &Node) -> bool {
    let (depth, width) = measure(a);
    depth > width && depth > b.depth
}

fn wider(a: &Node) -> u32 {
    let (depth, width) = measure(a);
    width.saturating_sub(depth)
}

// Flagged, with the evidence: both callers bind `(depth, width)` and the
// tail returns them in that order, but the early return spells
// `(width, depth)`. One of the two orders is a bug and the type accepts both.
fn clamped(n: &Node, max: u32) -> (u32, u32) {
    let depth = n.depth.min(max);
    let width = n.width.min(max);
    if depth == max {
        return (width, depth);
    }
    (depth, width)
}

fn clamped_area(n: &Node) -> u32 {
    let (depth, width) = clamped(n, 10);
    depth * width
}

fn clamped_ratio(n: &Node) -> u32 {
    let (depth, width) = clamped(n, 100);
    depth / width.max(1)
}

// Flagged, without the evidence: the only `(width, depth)` is a closure's
// return, which is the closure's value and not this function's.
fn folded(n: &Node) -> (u32, u32) {
    let clamp = |width: u32, depth: u32| {
        if width > depth {
            return (width, depth);
        }
        (depth, depth)
    };
    let depth = clamp(n.width, n.depth).1;
    (depth, n.width)
}

fn folded_area(n: &Node) -> u32 {
    let (depth, width) = folded(n);
    depth * width
}

fn folded_ratio(n: &Node) -> u32 {
    let (depth, width) = folded(n);
    depth / width.max(1)
}

// Flagged: a destructuring assignment names the members by what it assigns
// them to, `(depth, width)` like the `let` in the other caller.
fn grown(n: &Node) -> (u32, u32) {
    (n.depth + 2, n.width + 2)
}

fn grow(n: &mut Node) {
    (n.depth, n.width) = grown(n);
}

fn grown_area(n: &Node) -> u32 {
    let (depth, width) = grown(n);
    depth * width
}

// Flagged: the `Option` around the tuple is read through `?`, `if let`,
// `match` and `.unwrap()`, and every pattern that reaches the tuple names it
// `(key, value)`; `None` arms and `is_none()` never reach it.
fn split(s: &str) -> Option<(&str, &str)> {
    s.split_once('=')
}

fn key_len(s: &str) -> Option<usize> {
    if split(s).is_none() {
        return None;
    }
    let (key, value) = split(s)?;
    Some(key.len() + value.len())
}

fn has_value(s: &str) -> bool {
    if let Some((key, value)) = split(s) {
        return !key.is_empty() && !value.is_empty();
    }
    match split(s) {
        Some((key, value)) => key < value,
        None => false,
    }
}

fn value_of(s: &str) -> &str {
    let (key, value) = split(s).unwrap();
    if key.is_empty() { s } else { value }
}

struct Table {
    rows: Vec<(u8, u8)>,
}

struct Unbounded;

impl Table {
    // Flagged: an inherent method behind `Result`, read through `map_err`,
    // `?` and `let .. else`.
    fn bounds(&self) -> Result<(u8, u8), ()> {
        self.rows.first().copied().ok_or(())
    }

    fn span(&self) -> Result<u8, Unbounded> {
        let (low, high) = self.bounds().map_err(|()| Unbounded)?;
        Ok(high - low)
    }

    fn low(&self) -> u8 {
        let Ok((low, high)) = self.bounds() else {
            return 0;
        };
        low.min(high)
    }
}

// Fine: the members differ in type, so a transposed pattern or return is
// already a type error and the tuple loses nothing a struct would keep.
fn indexed(n: u32) -> (usize, u32) {
    (n as usize, n)
}

fn indexed_sum(n: u32) -> usize {
    let (index, value) = indexed(n);
    index + value as usize
}

fn indexed_gap(n: u32) -> usize {
    let (index, value) = indexed(n);
    index - value as usize
}

// Fine: the callers disagree (`(lo, hi)` here, `(start, end)` there), so
// neither pair is the members' name.
fn range(n: u32) -> (u32, u32) {
    (n, n + 10)
}

fn range_lo(n: u32) -> u32 {
    let (lo, hi) = range(n);
    lo.min(hi)
}

fn range_start(n: u32) -> u32 {
    let (start, end) = range(n);
    start.max(end)
}

// Fine: the callers disagree through destructuring assignments, `(near, far)`
// here and `(low, high)` there; the `lhs` rustc binds them to on the way is
// nobody's name for them.
fn assigned(n: u32) -> (u32, u32) {
    (n, n + 2)
}

fn assigned_near(n: u32) -> u32 {
    let (near, far);
    (near, far) = assigned(n);
    far - near
}

fn assigned_low(n: u32) -> u32 {
    let (low, high);
    (low, high) = assigned(n);
    low * high
}

// Fine: one caller keeps the tuple whole (`.0`), so the type is used as one.
fn halves(n: u32) -> (u32, u32) {
    (n / 2, n - n / 2)
}

fn halves_sum(n: u32) -> u32 {
    let (left, right) = halves(n);
    left + right
}

fn halves_left(n: u32) -> u32 {
    let (left, right) = halves(n);
    left * right + halves(n).0
}

// Fine: both destructuring sites sit in one function; one function's habit
// is not the crate's vocabulary.
fn corners(n: u32) -> (u32, u32) {
    (n, n * n)
}

fn corners_twice(n: u32) -> u32 {
    let (near, far) = corners(n);
    let (near2, far2) = {
        let (near, far) = corners(n + 1);
        (near, far)
    };
    near + far + near2 + far2
}

// Fine: a member is bound as `_` at one site, so that site names nothing.
fn parts(n: u32) -> (u32, u32) {
    (n, n + 1)
}

fn parts_head(n: u32) -> u32 {
    let (head, _) = parts(n);
    head
}

fn parts_both(n: u32) -> u32 {
    let (head, tail) = parts(n);
    head + tail
}

// Fine: the function is also passed as a value, so some destructuring sites
// are out of sight.
fn pair(n: u32) -> (u32, u32) {
    (n, n)
}

fn pair_a(n: u32) -> u32 {
    let (first, second) = pair(n);
    first + second
}

fn pair_b(n: u32) -> u32 {
    let (first, second) = pair(n);
    let f: fn(u32) -> (u32, u32) = pair;
    first * second + f(n).1
}

// Fine: a trait fixes the return type; the impl cannot change it alone.
trait Cursor {
    fn position(&self) -> (usize, usize);
}

impl Cursor for Table {
    fn position(&self) -> (usize, usize) {
        (self.rows.len(), 0)
    }
}

fn line_of(t: &Table) -> usize {
    let (line, column) = t.position();
    line + column
}

fn column_of(t: &Table) -> usize {
    let (line, column) = t.position();
    column.saturating_sub(line)
}

// Fine: exported, so callers outside the crate are invisible.
pub mod api {
    pub fn extent(n: usize) -> (usize, usize) {
        (n, n * 2)
    }
}

fn extent_a(n: usize) -> usize {
    let (offset, len) = api::extent(n);
    offset + len
}

fn extent_b(n: usize) -> usize {
    let (offset, len) = api::extent(n);
    offset * len
}

fn main() {
    let mut n = Node { depth: 1, width: 2 };
    grow(&mut n);
    let t = Table { rows: vec![(1, 2)] };
    let _ = (taller(&n, &n), wider(&n), clamped_area(&n), clamped_ratio(&n));
    let _ = (folded_area(&n), folded_ratio(&n), grown_area(&n));
    let _ = (key_len("a=b"), has_value("a=b"), value_of("a=b"));
    let _ = (t.span(), t.low());
    let _ = range_lo(1) + range_start(2) + halves_sum(3) + halves_left(4) + indexed_sum(5) as u32;
    let _ = assigned_near(3) + assigned_low(4);
    let _ = indexed_gap(9);
    let _ = corners_twice(5) + parts_head(6) + parts_both(7) + pair_a(8) + pair_b(9);
    let _ = line_of(&t) + column_of(&t) + extent_a(1) + extent_b(2);
}
