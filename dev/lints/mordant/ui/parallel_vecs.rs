// Sequence fields that only change length together and are read at one index are one Vec of a record.
#![allow(dead_code)]

use std::collections::VecDeque;

// Flagged: `names` and `ages` are pushed together, truncated together, and read at one index.
struct People {
    names: Vec<String>,
    ages: Vec<u32>,
    seen: usize,
}

impl People {
    // Both start empty, whatever the constructor is called.
    fn new() -> Self {
        Self {
            names: Vec::new(),
            ages: Vec::with_capacity(4),
            seen: 0,
        }
    }

    // `&mut` methods that cannot change a length are not writes.
    fn tidy(&mut self) {
        self.ages.reserve(8);
        self.names.sort();
        self.names.as_mut_slice().reverse();
    }

    fn add(&mut self, name: String, age: u32) {
        self.names.push(name);
        self.ages.push(age);
        self.seen += 1;
    }

    fn forget(&mut self, keep: usize) {
        self.names.truncate(keep);
        self.ages.truncate(keep);
    }

    fn describe(&self, i: usize) -> String {
        format!("{} is {}", self.names[i], self.ages[i])
    }

    // A `for` over `&mut self.ages` changes no length.
    fn birthday(&mut self) {
        for age in &mut self.ages {
            *age += 1;
        }
    }
}

// Flagged: three columns grown through a local binding and read by `zip`.
struct Table {
    keys: Vec<u32>,
    values: Vec<String>,
    phases: VecDeque<u8>,
}

fn insert(t: &mut Table, key: u32, value: String, phase: u8) {
    t.keys.push(key);
    t.values.push(value);
    t.phases.push_back(phase);
}

fn dump(t: &Table) {
    for (k, v) in t.keys.iter().zip(&t.values) {
        println!("{k} {v}");
    }
}

// Flagged: slices from one mark, through a field of `self`.
struct Tape {
    items: Vec<u8>,
    locs: Vec<u32>,
}

struct Builder {
    tape: Tape,
}

impl Builder {
    fn item(&mut self, item: u8, loc: u32) {
        self.tape.items.push(item);
        self.tape.locs.push(loc);
    }

    fn close(&mut self, mark: usize) -> usize {
        let n = self.tape.items[mark..].len() + self.tape.locs[mark..].len();
        self.tape.items.truncate(mark);
        self.tape.locs.truncate(mark);
        n
    }
}

// Fine: `errors` is also pushed alone, so it is not in step with `lines`.
struct Report {
    lines: Vec<String>,
    errors: Vec<String>,
}

impl Report {
    fn line(&mut self, l: String) {
        self.lines.push(l.clone());
        self.errors.push(l);
    }

    fn fail(&mut self, e: String) {
        self.errors.push(e);
    }

    fn pair(&self, i: usize) -> (&str, &str) {
        (&self.lines[i], &self.errors[i])
    }
}

// Fine: grown together but never read in step; each is consumed whole.
struct Rules {
    ltr: Vec<u8>,
    rtl: Vec<u8>,
}

impl Rules {
    fn add(&mut self, l: u8, r: u8) {
        self.ltr.push(l);
        self.rtl.push(r);
    }

    fn total(&self) -> usize {
        self.ltr.iter().map(|b| *b as usize).sum::<usize>() + self.rtl.len()
    }
}

// Fine: pushes to two different values are not a pair.
struct Lanes {
    xs: Vec<u8>,
    ys: Vec<u8>,
}

fn cross(a: &mut Lanes, b: &mut Lanes) {
    a.xs.push(1);
    b.ys.push(2);
}

fn lane(l: &Lanes, i: usize) -> u8 {
    l.xs[i] + l.ys[i]
}

// Fine: a `&mut` borrow handed to a function is a lone length change.
struct Swap {
    a: Vec<u8>,
    b: Vec<u8>,
}

impl Swap {
    fn both(&mut self, x: u8) {
        self.a.push(x);
        self.b.push(x);
    }

    fn reset(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.a)
    }

    fn sum(&self, i: usize) -> u8 {
        self.a[i] + self.b[i]
    }
}

// Fine: the two pushes sit in different match arms.
enum Side {
    Left,
    Right,
}

struct Split {
    left: Vec<u8>,
    right: Vec<u8>,
}

impl Split {
    fn put(&mut self, side: Side, x: u8) {
        match side {
            Side::Left => self.left.push(x),
            Side::Right => self.right.push(x),
        }
    }

    fn at(&self, i: usize) -> u8 {
        self.left[i] + self.right[i]
    }
}

// Fine: public fields of an exported struct can be grown by any crate.
pub struct Open {
    pub firsts: Vec<u8>,
    pub seconds: Vec<u8>,
}

pub fn open(o: &mut Open, i: usize) -> u8 {
    o.firsts.push(1);
    o.seconds.push(2);
    o.firsts[i] + o.seconds[i]
}

// Fine: `bytes` also grows alone through `io::Write`, a trait method the lint has no name for.
struct Log {
    bytes: Vec<u8>,
    marks: Vec<u8>,
}

impl Log {
    fn mark(&mut self, b: u8) {
        self.bytes.push(b);
        self.marks.push(b);
    }

    fn raw(&mut self, data: &[u8]) {
        use std::io::Write;
        self.bytes.write_all(data).unwrap();
    }

    fn at(&self, i: usize) -> u8 {
        self.bytes[i] + self.marks[i]
    }
}

// Fine: `text` is a `String` grown alone by `push_str`.
struct Masked {
    text: String,
    mask: Vec<bool>,
}

impl Masked {
    fn put(&mut self, c: char, m: bool) {
        self.text.push(c);
        self.mask.push(m);
    }

    fn word(&mut self, w: &str) {
        self.text.push_str(w);
    }

    fn at(&self, i: usize) -> bool {
        self.text.get(i..).is_some() && self.mask.get(i..).is_some()
    }
}

// Fine: `xs` shrinks alone through `pop_if`.
struct Popped {
    xs: Vec<u8>,
    ys: Vec<u8>,
}

impl Popped {
    fn add(&mut self, x: u8) {
        self.xs.push(x);
        self.ys.push(x);
    }

    fn trim(&mut self) {
        self.xs.pop_if(|x| *x == 0);
    }

    fn at(&self, i: usize) -> u8 {
        self.xs[i] + self.ys[i]
    }
}

// Fine: a custom sequence type grown alone through its own method.
struct List<T>(Vec<T>);

impl<T: Copy> List<T> {
    fn push(&mut self, t: T) {
        self.0.push(t);
    }

    fn append_slice(&mut self, s: &[T]) {
        self.0.extend_from_slice(s);
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, i: usize) -> Option<&T> {
        self.0.get(i)
    }
}

struct Lists {
    a: List<u8>,
    b: List<u8>,
}

impl Lists {
    fn add(&mut self, x: u8) {
        self.a.push(x);
        self.b.push(x);
    }

    fn bulk(&mut self, s: &[u8]) {
        self.a.append_slice(s);
    }

    fn at(&self, i: usize) -> u8 {
        *self.a.get(i).unwrap() + *self.b.get(i).unwrap()
    }
}

// Fine: `names` is cleared alone through a `&mut` binding split off `self`.
struct Destr {
    names: Vec<u8>,
    ages: Vec<u8>,
}

impl Destr {
    fn add(&mut self, x: u8) {
        self.names.push(x);
        self.ages.push(x);
    }

    fn only_names(&mut self) {
        let Self { names, .. } = self;
        names.clear();
    }

    fn at(&self, i: usize) -> u8 {
        self.names[i] + self.ages[i]
    }
}

// Fine: `xs` grows alone through a raw pointer to it.
struct Raw {
    xs: Vec<u8>,
    ys: Vec<u8>,
}

impl Raw {
    fn add(&mut self, x: u8) {
        self.xs.push(x);
        self.ys.push(x);
    }

    fn poke(&mut self) {
        let p = &raw mut self.xs;
        unsafe { (*p).push(0) };
    }

    fn at(&self, i: usize) -> u8 {
        self.xs[i] + self.ys[i]
    }
}

// Fine: born with `a` already longer than `b`.
struct Birth {
    a: Vec<u8>,
    b: Vec<u8>,
}

impl Birth {
    fn new(seed: Vec<u8>) -> Self {
        Birth {
            a: seed,
            b: Vec::new(),
        }
    }

    fn add(&mut self, x: u8) {
        self.a.push(x);
        self.b.push(x);
    }

    fn at(&self, i: usize) -> u8 {
        self.a[i] + self.b[i]
    }
}

fn main() {}
