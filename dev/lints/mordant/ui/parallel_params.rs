// Parameters several functions declare alike and hand on together are one value the types never name.

// Flagged: `w` and `h` are declared by `decode`, `scale` and `checksum`, and
// `decode` hands both unchanged to `scale`.
fn decode(src: &[u8], w: u32, h: u32) -> Vec<u8> {
    let out = scale(src, w, h, 2);
    let _ = checksum(w, h, &out);
    out
}

fn scale(src: &[u8], w: u32, h: u32, factor: u32) -> Vec<u8> {
    src.iter()
        .cycle()
        .take((w * h * factor) as usize)
        .copied()
        .collect()
}

fn checksum(w: u32, h: u32, px: &[u8]) -> u64 {
    u64::from(w) * u64::from(h) + px.len() as u64
}

// Flagged: `host`, `port` and `tls` through methods; the receiver is not one
// of them, and neither is the `&mut` sink, a context threaded through by
// design rather than part of the value.
#[derive(Clone, Copy)]
enum Tls {
    Off,
    On,
}

struct Conn {
    retries: u8,
}

impl Conn {
    fn send(&self, host: &str, port: u16, tls: Tls, sink: &mut Vec<String>) {
        for _ in 0..self.retries {
            self.dial(host, port, tls, &mut *sink);
        }
    }

    fn dial(&self, host: &str, port: u16, tls: Tls, sink: &mut Vec<String>) {
        log_target(host, port, tls, sink);
    }
}

fn log_target(host: &str, port: u16, tls: Tls, sink: &mut Vec<String>) {
    sink.push(format!("{host}:{port} {}", tls as u8));
}

// Fine: a builder and the log it reports to are both contexts, not data, so
// there is no value here to name however many functions pass them along.
struct Builder {
    depth: u32,
}

fn lower_block(builder: &mut Builder, log: &mut Vec<String>, stmts: &[u8]) {
    for s in stmts {
        lower_stmt(builder, log, *s);
    }
}

fn lower_stmt(builder: &mut Builder, log: &mut Vec<String>, stmt: u8) {
    builder.depth += 1;
    lower_expr(builder, log, stmt / 2);
    builder.depth -= 1;
}

fn lower_expr(builder: &mut Builder, log: &mut Vec<String>, expr: u8) {
    log.push(format!("{} at depth {}", expr, builder.depth));
}

// Fine: only `derive` and `mix` declare `key` and `salt`; the threshold is
// three functions.
fn derive(key: &[u8], salt: &[u8]) -> u8 {
    mix(key, salt) ^ 0x5c
}

fn mix(key: &[u8], salt: &[u8]) -> u8 {
    key.iter().chain(salt).fold(0, |a, b| a ^ b)
}

// Fine: three functions declare `lo` and `hi`, but only `clamp` hands the
// pair on, to `within`; `width` declaring the same names is not travel, so
// two functions are linked and the threshold is three.
fn clamp(v: i64, lo: i64, hi: i64) -> i64 {
    if within(v, lo, hi) { v } else { lo }
}

fn within(v: i64, lo: i64, hi: i64) -> bool {
    lo <= v && v <= hi
}

fn width(lo: i64, hi: i64) -> i64 {
    hi - lo
}

// Fine: `slice` renames what it receives, so `copy_range`'s call is a
// translation, not a hand-off; the other two never forward.
fn copy_range(buf: &[u8], from: usize, to: usize) -> Vec<u8> {
    slice(buf, from, to).to_vec()
}

fn slice(buf: &[u8], start: usize, end: usize) -> &[u8] {
    &buf[start..end]
}

fn check_range(buf: &[u8], from: usize, to: usize) -> bool {
    from <= to && to <= buf.len()
}

fn print_range(from: usize, to: usize) {
    println!("{from}..{to}");
}

// Fine: `resize` is a trait method and its impl, whose signatures the trait
// dictates; that leaves two functions declaring `rows` and `cols`.
trait Shape {
    fn resize(&mut self, rows: usize, cols: usize);
}

struct Grid {
    cells: Vec<u8>,
}

impl Shape for Grid {
    fn resize(&mut self, rows: usize, cols: usize) {
        self.cells = grid_cells(rows, cols);
    }
}

fn grid_cells(rows: usize, cols: usize) -> Vec<u8> {
    grid_fill(rows, cols, 0)
}

fn grid_fill(rows: usize, cols: usize, byte: u8) -> Vec<u8> {
    vec![byte; rows * cols]
}

// Fine: `invert` is also passed as a fn pointer, which pins its signature;
// without it only `blend` and `luma` declare `r`, `g`, `b`.
fn apply(px: &mut [u8; 3], op: fn(u8, u8, u8) -> u8) {
    px[0] = op(px[0], px[1], px[2]);
}

fn blend(r: u8, g: u8, b: u8) -> u8 {
    luma(r, g, b) / 2 + invert(r, g, b) / 2
}

fn luma(r: u8, g: u8, b: u8) -> u8 {
    r / 3 + g / 3 + b / 3
}

fn invert(r: u8, g: u8, b: u8) -> u8 {
    255 - luma(r, g, b)
}

// Fine: exported functions' signatures belong to callers this crate cannot
// see.
pub fn open(path: &str, mode: u32) -> usize {
    create(path, mode) + probe(path, mode)
}

pub fn create(path: &str, mode: u32) -> usize {
    path.len() + mode as usize
}

pub fn probe(path: &str, mode: u32) -> usize {
    path.len() ^ mode as usize
}

fn main() {
    let _ = decode(&[1, 2, 3], 4, 4);
    let mut log = Vec::new();
    Conn { retries: 2 }.send("h", 80, Tls::Off, &mut log);
    Conn { retries: 2 }.send("h", 443, Tls::On, &mut log);
    lower_block(&mut Builder { depth: 0 }, &mut log, b"xy");
    let _ = derive(b"k", b"s");
    let _ = (clamp(5, 0, 9), within(5, 0, 9), width(0, 9));
    let _ = (copy_range(b"abc", 0, 2), check_range(b"abc", 0, 2));
    print_range(0, 2);
    let mut g = Grid { cells: Vec::new() };
    g.resize(2, 2);
    let mut px = [1, 2, 3];
    apply(&mut px, invert);
    let _ = blend(1, 2, 3);
    let _ = open("/", 0o644);
}
