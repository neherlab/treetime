// One case per function. Each comment says the shape and whether it is reported, with the reason.
// "The part" is the run of statements the lint would move into a separate non-generic function.

// Reported: only `as_ref` and the drop of `bytes` use `B`. No help to take `&[u8]` instead:
// `checksum` is `pub`, so its signature is not this crate's to change.
pub fn checksum<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let src = bytes.as_ref();
    let mut acc = 17u32;
    let mut run = 0u32;
    for byte in src {
        let v = u32::from(*byte);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        if v & 1 == 0 {
            run += 1;
        } else {
            run = 0;
        }
        acc ^= run << 3;
    }
    (acc ^ (src.len() as u32)).rotate_left(7)
}

pub struct WithHeader<T> {
    pub payload: T,
    pub header: [u8; 4],
    pub declared_len: u32,
}

impl<T> WithHeader<T> {
    // Reported: `T` is the struct's parameter, but the method is still compiled once per `T`.
    pub fn header_word(&self) -> u32 {
        let header = self.header;
        let declared = u32::from_be(self.declared_len);
        let mut word = 0u32;
        let mut shift = 0u32;
        for b in header {
            word |= u32::from(b) << shift;
            shift += 8;
        }
        let padded = (declared + 3) & !3;
        if word > padded {
            word - padded
        } else {
            padded.wrapping_sub(word).wrapping_mul(3)
        }
    }
}

// Reported: a const parameter counts too. Past `as_slice` the body reads only a `&[u8]`.
pub fn fold_block<const N: usize>(block: [u8; N]) -> u32 {
    let bytes = block.as_slice();
    let mut lo = 1u32;
    let mut hi = 0u32;
    let mut i = 0usize;
    while i < bytes.len() {
        lo = (lo + u32::from(bytes[i])) % 65521;
        hi = (hi + lo) % 65521;
        i += 1;
    }
    let folded = (hi << 16) | lo;
    if folded % 2 == 0 { folded / 2 } else { folded.wrapping_mul(3) + 1 }
}

// Reported: reached only through `relay`, which is compiled twice. `#[inline]` does not stop it.
#[inline]
pub fn digits<S: AsRef<str>>(text: S) -> u32 {
    let s = text.as_ref();
    let mut total = 0u32;
    let mut weight = 1u32;
    let mut groups = 0u32;
    for ch in s.bytes() {
        if ch.is_ascii_digit() {
            total = total.wrapping_add(u32::from(ch - b'0') * weight);
            weight = weight.wrapping_mul(10);
        } else {
            if weight != 1 {
                groups += 1;
            }
            weight = 1;
        }
    }
    (total ^ weight.rotate_right(5)).wrapping_add(groups)
}

// Not reported: two statements per copy. This is the shape the lint asks for.
pub fn relay<S: AsRef<str>>(text: S) -> u32 {
    digits(text).wrapping_add(1)
}

// Reported, with the extra help: private, and uses `S` only to get a `&str`, so it could take `&str`.
fn digits_priv<S: AsRef<str>>(text: S) -> u32 {
    let s = text.as_ref();
    let mut vowels = 0u32;
    let mut longest = 0u32;
    let mut current = 0u32;
    for ch in s.bytes() {
        if matches!(ch, b'a' | b'e' | b'i' | b'o' | b'u') {
            vowels += 1;
            current += 1;
            if current > longest {
                longest = current;
            }
        } else {
            current = 0;
        }
    }
    (vowels << 8) | (longest & 0xff)
}

pub trait HasName {
    fn name(&self) -> &str;
}

pub struct Ann;
pub struct Bob;
pub struct Cy;

impl HasName for Ann {
    fn name(&self) -> &str {
        "ann"
    }
}

impl HasName for Bob {
    fn name(&self) -> &str {
        "Bob"
    }
}

impl HasName for Cy {
    fn name(&self) -> &str {
        "cy"
    }
}

// Reported. `tally` uses `X` only to get a `&str`, but it is `pub`, so no help to take `&str`.
pub fn tally<X: HasName>(x: X) -> u32 {
    let s = x.name();
    let bytes = s.as_bytes();
    let mut h = 2166136261u32;
    let mut upper = 0u32;
    let mut i = 0usize;
    while i < bytes.len() {
        let b = bytes[i];
        h = (h ^ u32::from(b)).wrapping_mul(16777619);
        if b.is_ascii_uppercase() {
            upper += 1;
        }
        i += 1;
    }
    (h ^ upper).wrapping_add(bytes.len() as u32)
}

// Reported: `if WIDE` splits the body. The first loop is reported, the note counts the second.
pub fn wide_sum<const WIDE: bool>(src: &[u8]) -> u32 {
    let mut acc = 1u32;
    let mut i = 0usize;
    while i < src.len() {
        let v = u32::from(src[i]);
        acc = acc.rotate_left(5) ^ v;
        acc = acc.wrapping_mul(3).wrapping_add(v >> 1);
        if acc & 0x10 == 0 {
            acc = acc.wrapping_add(i as u32);
        }
        i += 1;
    }
    if WIDE {
        acc = acc.swap_bytes();
    }
    let mut folded = acc;
    let mut j = 0usize;
    while j < src.len() {
        folded = folded.wrapping_add(u32::from(src[j]) << (j % 3));
        j += 2;
    }
    (folded ^ (src.len() as u32)).wrapping_add(7)
}

// Not reported: the long `else` branch runs only in the `FLAG = false` copy, so it is not shared.
pub fn arm_sum<const FLAG: bool>(src: &[u8]) -> u32 {
    let mut acc = 3u32;
    if FLAG {
        acc ^= src.len() as u32;
    } else {
        let mut i = 0usize;
        while i < src.len() {
            let v = u32::from(src[i]);
            acc = acc.rotate_left(5) ^ v;
            acc = acc.wrapping_mul(3).wrapping_add(v >> 1);
            if acc & 0x10 == 0 {
                acc = acc.wrapping_add(i as u32);
            }
            acc ^= acc >> 7;
            i += 1;
        }
        acc = acc.swap_bytes();
    }
    acc.wrapping_add(7)
}

#[derive(Clone, Copy)]
pub enum Coding {
    Plain,
    Hex,
    Sum,
    Xor,
}

pub const fn coding_of(tag: u8) -> Coding {
    match tag {
        0 => Coding::Plain,
        1 => Coding::Hex,
        2 => Coding::Sum,
        _ => Coding::Xor,
    }
}

// Not reported: the `match` is on a value computed from `TAG`, so each copy keeps one branch only.
pub fn encode_arm<const TAG: u8>(src: &[u8]) -> u32 {
    let coding = coding_of(TAG);
    match coding {
        Coding::Plain => {
            let mut acc = 0u32;
            for &b in src {
                acc = acc.wrapping_mul(31).wrapping_add(u32::from(b));
            }
            acc
        }
        Coding::Hex => {
            let mut acc = 1u32;
            for &b in src {
                acc = acc.rotate_left(4) ^ u32::from(b >> 4);
                acc = acc.rotate_left(4) ^ u32::from(b & 15);
            }
            acc
        }
        Coding::Sum => src.iter().map(|&b| u32::from(b)).sum::<u32>().wrapping_mul(7),
        Coding::Xor => {
            let mut acc = 0xffu32;
            for &b in src {
                acc ^= u32::from(b);
                acc = acc.rotate_right(1);
            }
            acc
        }
    }
}

pub const fn width_of(tag: u8) -> usize {
    match tag {
        0 => 1,
        1 => 2,
        2 => 4,
        _ => 8,
    }
}

// Not reported: the loop reads `width`, computed from `TAG`, which would stop being a constant.
pub fn encode_as<const TAG: u8>(src: &[u8]) -> u32 {
    let width = width_of(TAG);
    let mut acc = 0u32;
    let mut n = 0u32;
    for chunk in src.chunks(width) {
        let mut word = (width as u32).wrapping_mul(0x9e37_79b9);
        for &b in chunk {
            word = word.rotate_left(width as u32 + 3) ^ u32::from(b);
        }
        acc = acc.rotate_left(width as u32).wrapping_add(word);
        acc ^= acc >> (17 - width as u32);
        n = n.wrapping_add(width as u32);
    }
    acc ^ n
}

// Reported: `TAG` is stored into a byte of `buf`, which is run-time data, not a per-copy constant.
pub fn stamp_tag<const TAG: u8>(src: &[u8; 4]) -> [u8; 16] {
    let mut buf = [0u8; 16];
    buf[3] = TAG;
    buf[4] = src[0];
    buf[5] = src[1];
    buf[6] = src[2];
    buf[7] = src[3];
    buf[8] = buf[4] ^ buf[5];
    buf[9] = buf[6] ^ buf[7];
    buf[10] = buf[8].wrapping_add(buf[9]);
    buf[11] = buf[10].rotate_left(3);
    buf[12] = buf[11] ^ buf[4];
    buf[13] = buf[12].wrapping_mul(31);
    buf[14] = buf[13] ^ buf[5];
    buf[15] = buf[14].wrapping_add(buf[6]);
    buf[0] = buf[15];
    buf
}

// Reported, same part: `TAG` is only the index of a store, so `buf` is not a per-copy constant.
pub fn stamp_at<const TAG: u8>(src: &[u8; 4]) -> [u8; 16] {
    let mut buf = [0u8; 16];
    buf[TAG as usize] = 7;
    buf[4] = src[0];
    buf[5] = src[1];
    buf[6] = src[2];
    buf[7] = src[3];
    buf[8] = buf[4] ^ buf[5];
    buf[9] = buf[6] ^ buf[7];
    buf[10] = buf[8].wrapping_add(buf[9]);
    buf[11] = buf[10].rotate_left(3);
    buf[12] = buf[11] ^ buf[4];
    buf[13] = buf[12].wrapping_mul(31);
    buf[14] = buf[13] ^ buf[5];
    buf[15] = buf[14].wrapping_add(buf[6]);
    buf[0] = buf[15];
    buf
}

pub trait Sink {
    fn put(&mut self, byte: u8);
}

pub struct Count(pub u32);
pub struct Last(pub u8);

impl Sink for Count {
    fn put(&mut self, _byte: u8) {
        self.0 += 1;
    }
}

impl Sink for Last {
    fn put(&mut self, byte: u8) {
        self.0 = byte;
    }
}

// Not reported: the part is the loop body after `sink.put(b)`, so a call would run per iteration.
pub fn drain<S: Sink>(sink: &mut S, src: &[u8]) -> u32 {
    let mut acc = 5381u32;
    for &b in src {
        sink.put(b);
        let v = u32::from(b);
        acc = acc.rotate_left(5).wrapping_add(v);
        acc ^= v << 8;
        if v & 1 == 1 {
            acc = acc.wrapping_mul(33);
        } else {
            acc = acc.wrapping_sub(v >> 1);
        }
        acc = acc.rotate_right(v & 7) ^ 0x5bd1e995;
    }
    acc
}

// Reported: the `panic!` is inside the part. A block that panics is not a second exit.
pub fn strict_sum<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let src = bytes.as_ref();
    let mut acc = 1u32;
    let mut zeros = 0u32;
    for byte in src {
        let v = u32::from(*byte);
        if v == 0 {
            zeros += 1;
            if zeros > 3 {
                panic!("more than three zero bytes");
            }
        }
        acc = acc.wrapping_mul(37).wrapping_add(v);
        acc ^= acc >> 11;
    }
    (acc ^ zeros).rotate_right(3)
}

// Reported: one loop inside another. The outer loop is the larger part and is the one reported.
pub fn lattice<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let src = bytes.as_ref();
    let mut acc = 0x9e3779b9u32;
    let mut i = 0usize;
    while i < src.len() {
        let v = u32::from(src[i]);
        let mut bit = 0u32;
        while bit < 8 {
            if (v >> bit) & 1 == 1 {
                acc = acc.wrapping_mul(31) ^ bit;
            } else {
                acc = acc.rotate_left(3).wrapping_add(bit);
            }
            acc ^= acc >> 7;
            bit += 1;
        }
        acc ^= v;
        i += 1;
    }
    (acc ^ (src.len() as u32)).rotate_left(5)
}

// Not reported: `sink.put` runs every few statements, so no part between two calls is long enough.
pub fn render<S: Sink>(sink: &mut S, mut n: u32) -> u32 {
    n = n.wrapping_mul(2654435761);
    let hi = (n >> 24) as u8;
    sink.put(hi);
    n ^= n >> 13;
    let mid = (n >> 8) as u8;
    sink.put(mid);
    n = n.rotate_left(9).wrapping_add(40503);
    let lo = n as u8;
    sink.put(lo);
    n ^= u32::from(hi) | (u32::from(lo) << 16);
    sink.put((n % 251) as u8);
    n.count_ones() + u32::from(mid)
}

// Not reported: `f` is called in the middle of each iteration. The parts either side are too short.
pub fn scan<F: FnMut(u8)>(src: &[u8], mut f: F) -> u32 {
    let mut acc = 0u32;
    let mut prev = 0u8;
    for &b in src {
        let delta = b.wrapping_sub(prev);
        acc = acc.rotate_left(3) ^ u32::from(delta);
        f(delta);
        prev = b;
        acc = acc.wrapping_add(u32::from(b) << 2);
        if acc & 1 == 1 {
            acc ^= 0x9e3779b9;
        }
    }
    acc
}

// Not reported: `UP` is tested every few lines, so every part between the tests is too short.
pub fn stepped<const UP: bool>(mut n: u32) -> u32 {
    n = n.wrapping_mul(31).rotate_left(3);
    if UP {
        n = n.wrapping_add(7);
    }
    n ^= n >> 5;
    n = n.wrapping_mul(17);
    if UP {
        n = n.swap_bytes();
    }
    n = n.wrapping_sub(n >> 11) | 1;
    if UP {
        n ^= 0xa5a5;
    }
    n.wrapping_mul(3).rotate_right(2)
}

pub struct Wrap<T> {
    pub inner: T,
    pub len: u32,
    pub cap: u32,
}

impl<T> Wrap<T> {
    // Not reported: every statement reads a field through `&self`, and `&Wrap<T>` uses `T`.
    pub fn slack(&self) -> u32 {
        let mut n = self.cap.wrapping_sub(self.len);
        n = n.wrapping_mul(3) ^ self.cap;
        n = n.rotate_left(self.len & 15);
        n = n.wrapping_add(self.cap >> 2);
        n ^= self.len.count_ones();
        n = n.wrapping_mul(self.cap | 1);
        n = n.wrapping_sub(self.len >> 3);
        n.rotate_right(self.cap & 7)
    }
}

// Not reported: `#[inline(always)]` asks for a copy at every call site, so a shared call is out.
#[inline(always)]
pub fn mix<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let src = bytes.as_ref();
    let mut h = 0x811c9dc5u32;
    let mut i = 0usize;
    while i < src.len() {
        h = (h ^ u32::from(src[i])).wrapping_mul(0x01000193);
        h ^= h >> 15;
        h = h.rotate_left(1).wrapping_add(i as u32);
        i += 1;
    }
    h.wrapping_add(src.len() as u32).rotate_left(11)
}

// Not reported: called twice with `&str`, so there is one copy.
pub fn vowels<S: AsRef<str>>(text: S) -> u32 {
    let s = text.as_ref();
    let mut n = 0u32;
    let mut last_was = false;
    for ch in s.bytes() {
        let is = matches!(ch, b'a' | b'e' | b'i' | b'o' | b'u');
        if is && !last_was {
            n += 2;
        } else if is {
            n += 1;
        }
        last_was = is;
    }
    n.saturating_sub(1)
}

// Not reported: every statement moves, compares or copies a `T`.
pub fn largest<T: PartialOrd + Copy>(items: &[T], floor: T) -> T {
    let mut best = floor;
    for item in items {
        if *item > best {
            best = *item;
        }
    }
    best
}

// Not reported: two copies, but the shared part is under `generic-body-not-generic-min-statements`.
pub fn padded_len<B: AsRef<[u8]>>(bytes: B) -> usize {
    bytes.as_ref().len() * 2 + 1
}

// Not reported: already split. The generic part is one call and the work is in `spread_inner`.
pub fn spread<B: AsRef<[u8]>>(bytes: B) -> u32 {
    spread_inner(bytes.as_ref())
}

fn spread_inner(src: &[u8]) -> u32 {
    let mut min = u32::from(u8::MAX);
    let mut max = 0u32;
    for byte in src {
        let v = u32::from(*byte);
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    if src.is_empty() { 0 } else { (max - min) * 4 + max }
}

// Not reported: a lifetime parameter is erased before code generation, so there is one copy.
pub fn trailing_spaces<'a>(line: &'a str) -> &'a str {
    let bytes = line.as_bytes();
    let mut end = bytes.len();
    let mut seen = 0u32;
    while end > 0 && bytes[end - 1] == b' ' {
        end -= 1;
        seen += 1;
    }
    if seen > 4 { line } else { &line[..end] }
}

// Not reported: the arithmetic is in a closure, and closures are not measured.
pub fn weigh<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let step = |acc: u32, b: &u8| {
        let v = u32::from(*b);
        let mixed = acc.rotate_left(5) ^ v;
        let bumped = if v > 0x7f {
            mixed.wrapping_add(v * 3)
        } else {
            mixed.wrapping_sub(v + 11)
        };
        let folded = (bumped >> 16) ^ (bumped & 0xffff);
        if folded % 7 == 0 { folded | 1 } else { folded.wrapping_mul(40503) }
    };
    bytes.as_ref().iter().fold(5381, step)
}

// Reported as one part over all three loops, so the middle line does not split it. Its closures
// capture nothing, take and return plain integers, and their bodies do not name `S`.
pub fn map_in_part<S: Sink>(sink: &mut S, src: &[u8]) -> u32 {
    sink.put(0);
    let mut acc = 17u32;
    let mut run = 0u32;
    let mut odd = 0u32;
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        run = if v & 2 == 0 { run + 1 } else { 0 };
        odd ^= v.rotate_left(run & 7);
        acc ^= odd >> 3;
        acc = acc.wrapping_add(run);
    }
    let n = src.iter().map(|b| b ^ 0x55).fold(0u32, |n, v| n.wrapping_mul(33) ^ u32::from(v));
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(37).wrapping_add(v ^ n);
        run = if v & 4 == 0 { run + 2 } else { 1 };
        odd ^= v.rotate_left(run & 3);
        acc ^= odd >> 5;
        acc = acc.wrapping_add(run ^ odd);
    }
    acc ^ run ^ odd ^ n
}

// Reported as two parts, the loops before and after the middle line. Its first closure captures
// `key`, whose type uses `K`, so that line is compiled once per copy.
pub fn map_captures_generic<S: Sink, K: AsRef<[u8]>>(sink: &mut S, key: &K, src: &[u8]) -> u32 {
    sink.put(0);
    let mut acc = 17u32;
    let mut run = 0u32;
    let mut odd = 0u32;
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        run = if v & 2 == 0 { run + 1 } else { 0 };
        odd ^= v.rotate_left(run & 7);
        acc ^= odd >> 3;
        acc = acc.wrapping_add(run);
    }
    let n = src.iter().map(|b| b ^ key.as_ref()[0]).fold(0u32, |n, v| n.wrapping_mul(33) ^ u32::from(v));
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(37).wrapping_add(v ^ n);
        run = if v & 4 == 0 { run + 2 } else { 1 };
        odd ^= v.rotate_left(run & 3);
        acc ^= odd >> 5;
        acc = acc.wrapping_add(run ^ odd);
    }
    acc ^ run ^ odd ^ n
}

// Reported as two parts, like `map_captures_generic`. The closure captures nothing and takes and
// returns plain integers, but its body names `S`, so each copy differs.
pub fn map_body_names_param<S: Sink>(sink: &mut S, src: &[u8]) -> u32 {
    sink.put(0);
    let mut acc = 17u32;
    let mut run = 0u32;
    let mut odd = 0u32;
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        run = if v & 2 == 0 { run + 1 } else { 0 };
        odd ^= v.rotate_left(run & 7);
        acc ^= odd >> 3;
        acc = acc.wrapping_add(run);
    }
    let n = src.iter().map(|b| b ^ std::mem::size_of::<S>() as u8).fold(0u32, |n, v| n.wrapping_mul(33) ^ u32::from(v));
    for &b in src {
        let v = u32::from(b);
        acc = acc.wrapping_mul(37).wrapping_add(v ^ n);
        run = if v & 4 == 0 { run + 2 } else { 1 };
        odd ^= v.rotate_left(run & 3);
        acc ^= odd >> 5;
        acc = acc.wrapping_add(run ^ odd);
    }
    acc ^ run ^ odd ^ n
}

pub struct Sealed<T> {
    pub token: T,
    pub table: [u8; 8],
    pub seed: u32,
    pub word: u32,
}

// Reported: `deref` is reached through auto-deref and a deref coercion in `main`, one copy each.
impl<T> std::ops::Deref for Sealed<T> {
    type Target = u32;
    fn deref(&self) -> &u32 {
        let table = self.table;
        let mut acc = self.seed.rotate_left(1);
        let mut odd = 0u32;
        for b in table {
            let v = u32::from(b);
            acc = acc.rotate_left(3) ^ v;
            if v % 2 == 1 {
                odd += 1;
            }
        }
        let pick = (acc ^ odd) % 3 == 1 || odd > 4;
        if pick { &self.word } else { &self.seed }
    }
}

pub struct Record {
    pub name: Vec<u8>,
    pub count: u32,
    pub total: u32,
    pub sent: u32,
    pub flags: u32,
}

// Not reported: the part would take `rec: &mut Record` beside `name`, a borrow of `rec.name`.
pub fn recount<S: Sink>(sink: &mut S, rec: &mut Record, weight: u32) {
    let name: &[u8] = &rec.name;
    sink.put(rec.sent as u8);
    let len = name.len() as u32;
    let first = u32::from(name.first().copied().unwrap_or(0));
    let last = u32::from(name.last().copied().unwrap_or(0));
    rec.count = rec.count.wrapping_add(1);
    rec.total = rec.total.wrapping_add(len.wrapping_mul(weight));
    let mix = (first << 8 | last).rotate_left(rec.count & 31);
    rec.flags ^= mix;
    rec.flags = rec.flags.wrapping_mul(0x9e37_79b9) ^ weight;
    if rec.flags & 1 == 0 {
        rec.total = rec.total.rotate_left(3);
    } else {
        rec.count = rec.count.wrapping_add(len & 3);
    }
    rec.flags = rec.flags.wrapping_add(rec.total ^ rec.count);
    sink.put(last as u8);
    sink.put(name.len() as u8);
    rec.sent += 1;
}

// Not reported: `name = &rec.name` stays live across a part that would take `rec: &mut Record`.
pub fn restamp<S: Sink>(sink: &mut S, rec: &mut Record, weight: u32) {
    let name: &[u8] = &rec.name;
    let len = name.len() as u32;
    sink.put(len as u8);
    rec.count = rec.count.wrapping_add(1);
    rec.total = rec.total.wrapping_add(len.wrapping_mul(weight));
    let mix = (len << 8 | weight).rotate_left(rec.count & 31);
    rec.flags ^= mix;
    rec.flags = rec.flags.wrapping_mul(0x9e37_79b9) ^ weight;
    if rec.flags & 1 == 0 {
        rec.total = rec.total.rotate_left(3);
    } else {
        rec.count = rec.count.wrapping_add(len & 3);
    }
    rec.flags = rec.flags.wrapping_add(rec.total ^ rec.count);
    sink.put(name.first().copied().unwrap_or(0));
    rec.sent += 1;
}

// Reported: the tail from `rec.sent += 1` through the loop. The lines of `recount` come first.
// Checking their many refused candidates must not use up the limit before the tail is reached.
pub fn tail_after_refused_head<S: Sink>(sink: &mut S, rec: &mut Record, weight: u32, src: &[u8]) -> u32 {
    let name: &[u8] = &rec.name;
    sink.put(rec.sent as u8);
    let len = name.len() as u32;
    let first = u32::from(name.first().copied().unwrap_or(0));
    let last = u32::from(name.last().copied().unwrap_or(0));
    rec.count = rec.count.wrapping_add(1);
    rec.total = rec.total.wrapping_add(len.wrapping_mul(weight));
    let mix = (first << 8 | last).rotate_left(rec.count & 31);
    rec.flags ^= mix;
    rec.flags = rec.flags.wrapping_mul(0x9e37_79b9) ^ weight;
    if rec.flags & 1 == 0 {
        rec.total = rec.total.rotate_left(3);
    } else {
        rec.count = rec.count.wrapping_add(len & 3);
    }
    rec.flags = rec.flags.wrapping_add(rec.total ^ rec.count);
    sink.put(last as u8);
    sink.put(name.len() as u8);
    rec.sent += 1;
    let mut acc = 17u32;
    let mut run = 0u32;
    for byte in src {
        let v = u32::from(*byte);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        if v & 1 == 0 {
            run += 1;
        } else {
            run = 0;
        }
        acc ^= run << 3;
    }
    sink.put(acc as u8);
    (acc ^ (src.len() as u32)).rotate_left(7)
}

pub struct Cursor {
    pub data: Vec<u8>,
    pub pos: usize,
    pub reads: u32,
    pub sum: u32,
}

// Not reported: returning `chunk` would keep `*cur` borrowed while `cur.reads += 1` writes to it.
pub fn read_into<S: Sink>(sink: &mut S, cur: &mut Cursor, want: usize) {
    let start = cur.pos.min(cur.data.len());
    let end = cur.data.len().min(start + want);
    let chunk: &[u8] = &cur.data[start..end];
    let mut sum = cur.sum;
    for &b in chunk {
        sum = sum.rotate_left(5) ^ u32::from(b);
    }
    cur.sum = sum;
    cur.pos = end;
    let padded = end - start < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    sink.put(sum as u8);
    cur.reads += 1;
    sink.put(chunk.len() as u8);
}

pub trait Push {
    fn push(&mut self, bytes: &[u8]);
}

impl Push for Vec<u8> {
    fn push(&mut self, bytes: &[u8]) {
        self.extend_from_slice(bytes);
    }
}

impl Push for u32 {
    fn push(&mut self, bytes: &[u8]) {
        *self = self.wrapping_add(bytes.len() as u32);
    }
}

// Reported: only the part before `let chunk`, which produces integers, so nothing stays borrowed.
pub fn read_scanned<W: Push>(out: &mut W, cur: &mut Cursor, want: usize) {
    let start = cur.pos.min(cur.data.len());
    let avail = cur.data.len() - start;
    let take = avail.min(want).min(64);
    let end = start + take;
    let mut sum = cur.sum;
    let mut i = start;
    while i < end {
        sum = sum.rotate_left(5) ^ u32::from(cur.data[i]);
        i += 1;
    }
    cur.sum = sum;
    cur.pos = end;
    let chunk: &[u8] = &cur.data[start..end];
    let padded = take < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    cur.reads += 1;
    out.push(chunk);
    if padded {
        out.push(b"\0");
    }
}

// Not reported: as `read_into`, but `chunk` is pushed into `keep`, which then holds the borrow.
pub fn read_keep<'c, S: Sink>(sink: &mut S, cur: &'c mut Cursor, want: usize, keep: &mut Vec<&'c [u8]>) {
    let start = cur.pos.min(cur.data.len());
    let end = cur.data.len().min(start + want);
    let chunk: &[u8] = &cur.data[start..end];
    let mut sum = cur.sum;
    for &b in chunk {
        sum = sum.rotate_left(5) ^ u32::from(b);
    }
    cur.sum = sum;
    cur.pos = end;
    let padded = end - start < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    sink.put(sum as u8);
    keep.push(chunk);
    cur.reads += 1;
    sink.put(cur.reads as u8);
}

pub struct Kept<'k> {
    pub last: &'k [u8],
}

// Not reported: as `read_keep`, with the borrow stored through `out: &mut Kept` instead of pushed.
pub fn read_store<'c, S: Sink>(sink: &mut S, cur: &'c mut Cursor, want: usize, out: &mut Kept<'c>) {
    let start = cur.pos.min(cur.data.len());
    let end = cur.data.len().min(start + want);
    let chunk: &[u8] = &cur.data[start..end];
    let mut sum = cur.sum;
    for &b in chunk {
        sum = sum.rotate_left(5) ^ u32::from(b);
    }
    cur.sum = sum;
    cur.pos = end;
    let padded = end - start < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    sink.put(sum as u8);
    out.last = chunk;
    cur.reads += 1;
    sink.put(cur.reads as u8);
}

// Not reported: returning `chunk` keeps the `Box` borrowed while `cur.reads += 1` writes to it.
pub fn peek_boxed<S: Sink>(sink: &mut S, mut cur: Box<Cursor>, want: usize) -> Box<Cursor> {
    let start = cur.pos.min(cur.data.len());
    let end = cur.data.len().min(start + want);
    let chunk: &[u8] = &cur.data[start..end];
    let mut sum = cur.sum;
    for &b in chunk {
        sum = sum.rotate_left(5) ^ u32::from(b);
    }
    sink.put(sum as u8);
    cur.reads += 1;
    sink.put(chunk.len() as u8);
    cur
}

pub struct Src {
    pub buf: Vec<u8>,
    pub pos: usize,
}

pub struct Parser<'a> {
    pub src: &'a mut Src,
    pub depth: u32,
    pub sum: u32,
}

// Not reported: returning `tok` would keep all of `*p` mutably borrowed while `p.depth` is read.
pub fn next_token<S: Sink>(sink: &mut S, p: &mut Parser<'_>, want: usize) {
    let start = p.src.pos.min(p.src.buf.len());
    let end = p.src.buf.len().min(start + want);
    let tok: &[u8] = &p.src.buf[start..end];
    let mut sum = p.sum;
    for &b in tok {
        sum = sum.rotate_left(5) ^ u32::from(b);
    }
    p.src.pos = end;
    let padded = end - start < want;
    if padded {
        sum ^= 0xa5a5_a5a5;
    }
    sink.put(sum as u8);
    sink.put(p.depth as u8);
    sink.put(tok.len() as u8);
}

// Reported: ends at the `if`. Including `cur.reads += 1` would return its unfinished checked add.
pub fn refill<S: Sink>(sink: &mut S, cur: &mut Cursor, want: usize) {
    let start = cur.pos.min(cur.data.len());
    let avail = cur.data.len() - start;
    let take = avail.min(want).min(64);
    let end = start + take;
    let mut sum = cur.sum;
    let mut i = start;
    while i < end {
        sum = sum.rotate_left(5) ^ u32::from(cur.data[i]);
        i += 1;
    }
    cur.sum = sum;
    let padded = take < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    cur.pos = end;
    cur.reads += 1;
    sink.put((end - start) as u8);
    if padded {
        sink.put(0);
    }
}

// Reported, same part: two `+=` in a row, and both blocks are left out for the same reason.
pub fn refill_twice<S: Sink>(sink: &mut S, cur: &mut Cursor, want: usize) {
    let start = cur.pos.min(cur.data.len());
    let avail = cur.data.len() - start;
    let take = avail.min(want).min(64);
    let end = start + take;
    let mut sum = cur.sum;
    let mut i = start;
    while i < end {
        sum = sum.rotate_left(5) ^ u32::from(cur.data[i]);
        i += 1;
    }
    cur.sum = sum;
    cur.pos = end;
    let padded = take < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    cur.sum += 3;
    cur.reads += 1;
    sink.put((end - start) as u8);
    if padded {
        sink.put(0);
    }
}

// Reported, same part: the block of `let total = cur.reads + 1` is left out whole the same way.
pub fn refill_total<S: Sink>(sink: &mut S, cur: &mut Cursor, want: usize) -> u32 {
    let start = cur.pos.min(cur.data.len());
    let avail = cur.data.len() - start;
    let take = avail.min(want).min(64);
    let end = start + take;
    let mut sum = cur.sum;
    let mut i = start;
    while i < end {
        sum = sum.rotate_left(5) ^ u32::from(cur.data[i]);
        i += 1;
    }
    cur.sum = sum;
    cur.pos = end;
    let padded = take < want;
    if padded {
        cur.sum ^= 0xa5a5_a5a5;
    }
    let total = cur.reads + 1;
    sink.put((end - start) as u8);
    if padded {
        sink.put(0);
    }
    total
}

fn consume(text: String) -> u32 {
    text.len() as u32
}

// Not reported: `drop(guard)` runs on one path only, so the part would need a hidden drop flag.
pub fn settle<S: Sink>(sink: &mut S, lock: &std::sync::Mutex<u32>, seed: u32) -> u32 {
    let guard = lock.lock().unwrap();
    let base = *guard;
    sink.put(base as u8);
    let mut acc = seed ^ base;
    acc = acc.rotate_left(5).wrapping_mul(31);
    acc ^= acc >> 7;
    acc = acc.wrapping_add(base | 1);
    if acc & 1 == 0 {
        drop(guard);
    }
    acc = acc.wrapping_mul(2654435761);
    acc ^= acc >> 13;
    acc = acc.rotate_left(9).wrapping_add(seed);
    sink.put(acc as u8);
    acc
}

// Not reported: `text` is made on one path only, so the part would return a hidden drop flag.
pub fn label<S: Sink>(sink: &mut S, seed: u32) -> u32 {
    let text: String;
    sink.put(seed as u8);
    let mut acc = seed;
    acc = acc.rotate_left(5).wrapping_mul(31);
    acc ^= acc >> 7;
    acc = acc.wrapping_add(seed | 1);
    if acc & 1 == 0 {
        text = String::from("even");
        acc ^= text.len() as u32;
    }
    acc = acc.wrapping_mul(2654435761);
    acc ^= acc >> 13;
    acc = acc.rotate_left(9).wrapping_add(seed);
    sink.put(acc as u8);
    acc
}

// Reported: the part stops before the `return`, where the hidden drop flag of `text` is tested.
pub fn unlabel<S: Sink>(sink: &mut S, seed: u32, keep: bool) -> u32 {
    let text = String::from("kept");
    let mut acc = seed;
    if !keep {
        acc ^= consume(text);
    }
    sink.put(acc as u8);
    acc = acc.rotate_left(5).wrapping_mul(31);
    acc ^= acc >> 7;
    acc = acc.wrapping_add(seed | 1);
    acc = acc.rotate_right(3).wrapping_sub(40503);
    acc ^= acc >> 11;
    acc = acc.wrapping_mul(2654435761);
    acc ^= acc >> 13;
    acc = acc.rotate_left(9).wrapping_add(seed);
    acc = acc.wrapping_sub(acc >> 5) | 1;
    acc ^= acc >> 3;
    acc.wrapping_mul(16777619).rotate_right(seed & 7)
}

// Reported: the part takes `text` by value and may drop it. The hidden drop flag is not listed.
pub fn handoff<S: Sink>(sink: &mut S, seed: u32) -> u32 {
    let text = String::from("payload");
    sink.put(seed as u8);
    let mut acc = seed;
    acc = acc.rotate_left(5).wrapping_mul(31);
    acc ^= acc >> 7;
    acc = acc.wrapping_add(seed | 1);
    if acc & 1 == 0 {
        acc ^= consume(text);
    }
    acc = acc.wrapping_mul(2654435761);
    acc = acc.rotate_left(9).wrapping_add(seed);
    if acc & 2 == 0 {
        acc ^= 0x9e3779b9;
    }
    sink.put(acc as u8);
    acc
}

pub trait Emit {
    fn emit(&mut self, bytes: &[u8]);
}

impl Emit for Count {
    fn emit(&mut self, bytes: &[u8]) {
        self.0 += bytes.len() as u32;
    }
}

impl Emit for Last {
    fn emit(&mut self, bytes: &[u8]) {
        if let Some(&b) = bytes.last() {
            self.0 = b;
        }
    }
}

// Not reported: the part would return `view`, which borrows `buf`, a local made inside the part.
pub fn own_view<E: Emit>(out: &mut E, seed: u32) {
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    let view = &buf[..n];
    out.emit(view);
}

#[repr(C)]
pub struct Iov {
    pub base: *mut u8,
    pub len: usize,
}

pub struct Extent {
    pub first: usize,
    pub len: usize,
}

pub trait EmitRaw {
    fn emit_iov(&mut self, iov: &Iov, ready: usize);
    fn emit_non_null(&mut self, bytes: std::ptr::NonNull<[u8; 8]>, len: usize);
    fn emit_waker(&mut self, waker: std::task::RawWaker, len: usize);
    fn emit_extent(&mut self, extent: &Extent, ready: usize);
}

impl EmitRaw for Count {
    fn emit_iov(&mut self, iov: &Iov, ready: usize) {
        self.0 += (iov.len + ready) as u32;
    }
    fn emit_non_null(&mut self, _: std::ptr::NonNull<[u8; 8]>, len: usize) {
        self.0 += len as u32;
    }
    fn emit_waker(&mut self, _: std::task::RawWaker, len: usize) {
        self.0 += len as u32;
    }
    fn emit_extent(&mut self, extent: &Extent, ready: usize) {
        self.0 += (extent.first + extent.len + ready) as u32;
    }
}

impl EmitRaw for Last {
    fn emit_iov(&mut self, iov: &Iov, ready: usize) {
        self.0 = (iov.len + ready) as u8;
    }
    fn emit_non_null(&mut self, _: std::ptr::NonNull<[u8; 8]>, len: usize) {
        self.0 = len as u8;
    }
    fn emit_waker(&mut self, _: std::task::RawWaker, len: usize) {
        self.0 = len as u8;
    }
    fn emit_extent(&mut self, extent: &Extent, ready: usize) {
        self.0 = (extent.len + ready) as u8;
    }
}

// Not reported: the part would return `iov`, whose `base` field is the address of `buf`, a local made inside the part.
pub fn own_iov<E: EmitRaw>(out: &mut E, seed: u32) {
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let mut buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    let iov = Iov { base: buf.as_mut_ptr(), len: n };
    let ready = iov.len.min(4);
    out.emit_iov(&iov, ready);
}

// Not reported: the part would return `bytes`, a `NonNull` that is the address of `buf`, a local made inside the part.
pub fn own_non_null<E: EmitRaw>(out: &mut E, seed: u32) {
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let mut buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    let bytes = std::ptr::NonNull::from(&mut buf);
    out.emit_non_null(bytes, n);
}

static IDLE: std::task::RawWakerVTable =
    std::task::RawWakerVTable::new(|data| std::task::RawWaker::new(data, &IDLE), |_| {}, |_| {}, |_| {});

// Not reported: the part would return `waker`, which keeps the address of `buf` in a private field.
// `RawWaker` is from another crate but has no `Drop` impl, so it does not own what it points to.
pub fn own_waker<E: EmitRaw>(out: &mut E, seed: u32) {
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    let waker = std::task::RawWaker::new(buf.as_ptr().cast(), &IDLE);
    out.emit_waker(waker, n);
}

// Reported: the part borrows `buf` but returns `extent`, whose integer fields cannot hold that address.
pub fn own_extent<E: EmitRaw>(out: &mut E, seed: u32) {
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    let extent = Extent { first: usize::from(buf[0]), len: buf[..n].len() };
    let ready = extent.len.min(4);
    out.emit_extent(&extent, ready);
}

// Not reported: the borrow of `buf` is pushed into `parts`, which is read after the part.
pub fn parked_view<E: Emit>(out: &mut E, seed: u32) {
    let mut parts: Vec<&[u8]> = Vec::with_capacity(2);
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    parts.push(&buf[..n]);
    parts.push(b";");
    for part in &parts {
        out.emit(part);
    }
}

fn stash<'a>(slot: std::rc::Rc<std::cell::Cell<Option<&'a [u8]>>>, view: &'a [u8]) {
    slot.set(Some(view));
}

// Not reported: the borrow of `buf` is stored through an `Rc` whose `Cell` is read after the part.
pub fn rc_parked<E: Emit>(out: &mut E, seed: u32) {
    let slot = std::rc::Rc::new(std::cell::Cell::new(None));
    let slot2 = std::rc::Rc::clone(&slot);
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    stash(slot2, &buf[..n]);
    if let Some(v) = slot.get() {
        out.emit(v);
    }
}

pub struct Both<'a, 'b> {
    pub v: &'a [u8],
    pub p: &'b mut Vec<&'a [u8]>,
}

impl Both<'_, '_> {
    fn go(self) {
        self.p.push(self.v);
    }
}

// Not reported: `&buf[..n]` and `&mut parts` go into one `Both` value, refused like the call.
pub fn single_arg_park<E: Emit>(out: &mut E, seed: u32) {
    let mut parts: Vec<&[u8]> = Vec::with_capacity(2);
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    Both {
        v: &buf[..n],
        p: &mut parts,
    }
    .go();
    parts.push(b";");
    for part in &parts {
        out.emit(part);
    }
}

fn push_pair<'a>(pair: (&'a [u8], &mut Vec<&'a [u8]>)) {
    pair.1.push(pair.0);
}

// Not reported: the same through a tuple argument.
pub fn tuple_arg_park<E: Emit>(out: &mut E, seed: u32) {
    let mut parts: Vec<&[u8]> = Vec::with_capacity(2);
    let x = seed.wrapping_mul(0x9e37_79b9);
    let y = x.rotate_left(7) ^ seed;
    let z = y.swap_bytes().wrapping_add(x);
    let n = (z & 7) as usize;
    let buf = [
        x as u8,
        y as u8,
        z as u8,
        (x ^ y) as u8,
        (y & z) as u8,
        (x | z) as u8,
        !(x as u8),
        (x ^ y ^ z) as u8,
    ];
    push_pair((&buf[..n], &mut parts));
    parts.push(b";");
    for part in &parts {
        out.emit(part);
    }
}

fn advance(cur: &mut &mut [u8], bytes: &[u8]) {
    let (head, tail) = std::mem::take(cur).split_at_mut(bytes.len());
    head.copy_from_slice(bytes);
    *cur = tail;
}

// Reported: only from `start - len` on. Larger parts would return or alias a borrow of `buf`.
pub fn trim_stamp<E: Emit>(out: &mut E, secs: u32, frac: u32) -> usize {
    let mut buf = [b' '; 32];
    let mut cursor = &mut buf[4..];
    let start = cursor.len();
    advance(&mut cursor, &[b'0' | secs as u8, b'.', frac as u8, b'0']);
    let len = cursor.len();
    let n = start - len;
    let text = &buf[..4 + n];
    let mut end = text.len();
    while end > 0 && text[end - 1] == b' ' {
        end -= 1;
    }
    let mut begin = 0;
    while begin < end && text[begin] == b' ' {
        begin += 1;
    }
    if end - begin > 2 && text[end - 1] == b'0' && text[end - 2] == b'0' {
        end -= 2;
    }
    let shown = &text[begin..end];
    out.emit(shown);
    end - begin
}

pub struct FrameHeader {
    pub length: u32,
    pub kind: u8,
    pub flags: u8,
    pub stream: u32,
}

impl FrameHeader {
    // Reported: fills `buf` from `&self` and returns it by value, so nothing stays borrowed.
    pub fn write<E: Emit>(&self, out: &mut E) -> usize {
        let mut buf = [0u8; 9];
        buf[0] = (self.length >> 16) as u8;
        buf[1] = (self.length >> 8) as u8;
        buf[2] = self.length as u8;
        buf[3] = self.kind;
        buf[4] = self.flags;
        buf[5..9].copy_from_slice(&self.stream.to_be_bytes());
        out.emit(&buf);
        buf.len()
    }
}

pub struct Stats {
    pub dots: u32,
    pub longest: usize,
}

// Reported: only `let start` through `stats.dots = ..`. Larger parts would keep a borrow of `buf`.
pub fn print_stamp<E: Emit>(out: &mut E, stats: &mut Stats, secs: u64, frac: u32) {
    use std::io::Write as _;
    let mut buf = [0u8; 40];
    let mut cursor = &mut buf[..];
    let start = cursor.len();
    let days = secs / 86_400;
    let rem = secs - days * 86_400;
    let hours = (rem / 3_600) as u32;
    let minutes = ((rem / 60) % 60) as u32;
    let seconds = (rem % 60) as u32;
    let year_ish = 1970 + (days * 400 / 146_097) as u32;
    let day_ish = (days - u64::from(year_ish - 1970) * 365) as u32 % 366;
    let millis = frac / 1_000_000;
    let micros = (frac / 1_000) % 1_000;
    let packed = (hours << 12) | (minutes << 6) | seconds;
    stats.dots = stats.dots.wrapping_add(packed ^ day_ish ^ micros);
    let _ = write!(cursor, "{year_ish:04}-{day_ish:03}T{hours:02}:{minutes:02}:{seconds:02}");
    if frac != 0 {
        let _ = write!(cursor, ".{millis:03}{micros:03}");
    }
    let _ = write!(cursor, "Z");
    let n = start - cursor.len();
    let text = &buf[..n];
    stats.longest = stats.longest.max(n);
    out.emit(text);
}

pub struct Settings<M> {
    pub marker: M,
    pub width: u32,
    pub height: u32,
    pub depth: u32,
    pub gamma: u32,
    pub seed: u32,
    pub limit: u32,
    pub floor: u32,
    pub scale: u32,
}

// Reported: the eight field initializers. The note names the unnamed `u32`s they produce by field.
pub fn settings_for<M>(marker: M, base: u32) -> Settings<M> {
    Settings {
        width: base.wrapping_mul(3).rotate_left(2) ^ 0x51,
        height: base.wrapping_mul(5).rotate_left(3) ^ 0x52,
        depth: base.wrapping_mul(7).rotate_left(5) ^ 0x53,
        gamma: base.wrapping_mul(11).rotate_left(7) ^ 0x54,
        seed: base.wrapping_mul(13).rotate_left(11) ^ 0x55,
        limit: base.wrapping_mul(17).rotate_left(13) ^ 0x56,
        floor: base.wrapping_mul(19).rotate_left(17) ^ 0x57,
        scale: base.wrapping_mul(23).rotate_left(19) ^ 0x58,
        marker,
    }
}

pub struct Printer {
    pub col: u32,
    pub out: Vec<u8>,
    pub errors: u32,
}

impl Printer {
    fn fail(&mut self) -> u32 {
        self.errors += 1;
        self.errors
    }

    // Reported: `self.col = ..` through `self.out.len()`, producing the `usize` the `if` tests.
    pub fn write_str<S: AsRef<[u8]>>(&mut self, s: S) -> Result<(), u32> {
        let s = s.as_ref();
        self.col = self.col.wrapping_add(s.len() as u32);
        let mut h = self.col;
        for &b in s {
            h = (h ^ u32::from(b)).wrapping_mul(0x01000193);
        }
        self.out.extend_from_slice(s);
        self.out.push((h & 0x7f) as u8);
        if self.out.len() > 4096 {
            return Err(self.fail());
        }
        Ok(())
    }
}

pub struct Ledger<T> {
    pub tag: T,
    pub counts: [u32; 6],
    pub total: u32,
}

// Reported: nothing calls `drop` by name. The two `Ledger`s in `main` run it, one copy each.
impl<T> Drop for Ledger<T> {
    fn drop(&mut self) {
        let counts = self.counts;
        let floor = self.total.min(60);
        let mut sum = 0u32;
        let mut peak = 0u32;
        let mut at = 0u32;
        let mut i = 0u32;
        for c in counts {
            sum = sum.wrapping_add(c);
            if c > peak {
                peak = c;
                at = i;
            }
            i += 1;
        }
        let spread = peak.saturating_sub(sum / 6).wrapping_add(floor);
        let total = if spread > at { sum ^ spread } else { sum.wrapping_mul(at + 1) };
        self.total = total;
    }
}

// Not reported: the same body as `checksum`, but a macro wrote the function.
macro_rules! make_summer {
    ($name:ident) => {
        pub fn $name<B: AsRef<[u8]>>(bytes: B) -> u32 {
            let src = bytes.as_ref();
            let mut acc = 17u32;
            let mut run = 0u32;
            for byte in src {
                let v = u32::from(*byte);
                acc = acc.wrapping_mul(31).wrapping_add(v);
                if v & 1 == 0 {
                    run += 1;
                } else {
                    run = 0;
                }
                acc ^= run << 3;
            }
            (acc ^ (src.len() as u32)).rotate_left(7)
        }
    };
}

make_summer!(macro_checksum);

pub const VERBOSE: bool = true;

pub struct Scope {
    pub tag: &'static str,
    pub on: bool,
}

impl Scope {
    pub fn visible(&self) -> bool {
        self.on
    }
    pub fn log(&self, args: std::fmt::Arguments<'_>) {
        if self.on {
            eprintln!("{args}");
        }
    }
}

pub static REQUESTS: Scope = Scope { tag: "req", on: true };

// A logging macro. One use expands to about fifty statements that use no parameter.
macro_rules! trace {
    ($scope:path, $fmt:expr $(, $arg:expr)* $(,)?) => {
        if VERBOSE && $scope.visible() {
            if $scope.tag.len() > 2 {
                $scope.log(format_args!(concat!("\x1b[2m[{}]\x1b[0m ", $fmt, "{}"), $scope.tag, $($arg,)* "\n"));
            } else {
                $scope.log(format_args!(concat!("[{}] ", $fmt, "{}"), $scope.tag, $($arg,)* "\n"));
            }
        }
    };
}

pub struct Conn<const TLS: bool> {
    pub id: u32,
    pub open: bool,
}

impl<const TLS: bool> Conn<TLS> {
    // Not reported: the only long part is the `trace!` expansion. Macro statements do not count.
    pub fn close(&mut self) -> u32 {
        trace!(REQUESTS, "close");
        if !self.open {
            return self.id;
        }
        self.open = false;
        self.id.rotate_left(if TLS { 3 } else { 5 })
    }
}

// Not reported: both uses are `const` initializers, so no copy is compiled into the binary.
pub const fn const_fold<const N: usize>(block: [u8; N]) -> u32 {
    let bytes = block.as_slice();
    let mut lo = 1u32;
    let mut hi = 0u32;
    let mut i = 0usize;
    while i < bytes.len() {
        lo = (lo + bytes[i] as u32) % 65521;
        hi = (hi + lo) % 65521;
        i += 1;
    }
    let folded = (hi << 16) | lo;
    if folded % 2 == 0 { folded / 2 } else { folded.wrapping_mul(3) + 1 }
}

pub const FOLD_TWO: u32 = const_fold([1u8, 2]);
pub const FOLD_THREE: u32 = const_fold([1u8, 2, 3]);

fn evens(src: &[u8]) -> impl Iterator<Item = u8> + '_ {
    src.iter().copied().filter(|b| b & 1 == 0)
}

// Not reported: the part would take `it`, whose `impl Iterator` type no signature can name.
pub fn fold_evens<S: Sink>(sink: &mut S, src: &[u8]) -> u32 {
    let mut it = evens(src);
    sink.put(0);
    let mut acc = 17u32;
    let mut run = 0u32;
    while let Some(b) = it.next() {
        let v = u32::from(b);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        run = if v & 2 == 0 { run + 1 } else { 0 };
    }
    acc ^= run << 3;
    acc = (acc ^ run).rotate_left(7);
    sink.put(acc as u8);
    acc
}

pub struct Row {
    pub a: u32,
    pub b: u32,
    pub c: u32,
}

// Reported: the `else` arm uses `S`, so the part is the `if let` block. The note starts at that
// block, not at `Some(row)`: `row` is bound after `rows.get`, which the part only reads.
pub fn hash_row<S: Sink>(sink: &mut S, rows: &[Row], key: u32) -> u32 {
    let mut acc = key;
    if let Some(row) = rows.get(key as usize) {
        let v = row.a.rotate_left(3);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        acc ^= acc >> 7;
        acc = acc.wrapping_add(row.b | 1);
        acc = acc.rotate_left(5);
        acc ^= row.c.wrapping_mul(0x9e37_79b9);
        acc = acc.wrapping_sub(v >> 3);
        acc ^= acc >> 11;
        acc = acc.wrapping_mul(0x85eb_ca6b);
        acc ^= acc >> 13;
        acc = acc.wrapping_add(row.a ^ row.b);
        acc ^= acc >> 16;
        acc = acc.rotate_left(row.c & 31);
        acc ^= 0x5555_aaaa;
    } else {
        sink.put(0);
    }
    sink.put(acc as u8);
    acc
}

pub struct Dims {
    pub w: u32,
    pub h: u32,
}

impl Dims {
    pub fn split(&self) -> (u32, u32) {
        (self.w, self.h)
    }
}

pub struct Canvas<T> {
    pub tag: T,
    pub dims: Dims,
}

impl<T> Canvas<T> {
    // Reported: `self.dims.split()` reads `self`, a `&Canvas<T>`, so the part only reads its
    // result. The note starts at `w.rotate_left(3)`, not at `(w, h)`, which is bound after it.
    pub fn area<S: Sink>(&self, sink: &mut S) -> u32 {
        let (w, h) = self.dims.split();
        let mut acc = w.rotate_left(3);
        acc = acc.wrapping_mul(31).wrapping_add(h);
        acc ^= acc >> 7;
        acc = acc.wrapping_add(h | 1);
        acc = acc.rotate_left(5);
        acc ^= w.wrapping_mul(0x9e37_79b9);
        acc = acc.wrapping_sub(h >> 3);
        acc ^= acc >> 11;
        acc = acc.wrapping_mul(0x85eb_ca6b);
        acc ^= acc >> 13;
        acc = acc.wrapping_add(w ^ h);
        acc ^= acc >> 16;
        acc = acc.rotate_left(h & 31);
        acc ^= 0x5555_aaaa;
        sink.put(acc as u8);
        acc
    }
}

fn main() {
    let _ = checksum("abc");
    let _ = checksum(vec![1u8, 2, 3]);
    let _ = checksum([9u8; 4]);
    let small = WithHeader { payload: 1u8, header: [1, 2, 3, 4], declared_len: 6 };
    let wide = WithHeader { payload: "wide", header: [4, 3, 2, 1], declared_len: 9 };
    let _ = small.header_word() + wide.header_word();
    let _ = fold_block([1u8, 2]);
    let _ = fold_block([1u8, 2, 3]);
    let _ = relay("12a3");
    let _ = relay(String::from("45"));
    let _ = digits_priv("12a3") + digits_priv(String::from("45"));
    let _ = tally(Ann) + tally(Bob) + tally(Cy);
    let _ = wide_sum::<true>(b"wide") + wide_sum::<false>(b"narrow");
    let _ = arm_sum::<true>(b"this") + arm_sum::<false>(b"that");
    let _ = encode_arm::<0>(b"p") + encode_arm::<1>(b"h") + encode_arm::<2>(b"s");
    let _ = encode_arm::<3>(b"x");
    let _ = encode_as::<0>(b"ab") ^ encode_as::<1>(b"cd") ^ encode_as::<2>(b"ef");
    let _ = encode_as::<3>(b"gh");
    let _ = stamp_tag::<0>(b"abcd")[0] ^ stamp_tag::<1>(b"cdef")[0];
    let _ = stamp_at::<2>(b"abcd")[0] ^ stamp_at::<3>(b"cdef")[0];
    let mut count = Count(0);
    let mut last = Last(0);
    let _ = drain(&mut count, b"abc") + drain(&mut last, b"de");
    let _ = strict_sum("s") + strict_sum([1u8, 0]);
    let _ = lattice("l") + lattice(vec![3u8]);
    let _ = render(&mut count, 3) + render(&mut last, 4);
    let _ = scan(b"up", |b| count.0 += u32::from(b));
    let _ = scan(b"down", |b| last.0 = b);
    let _ = stepped::<true>(5) + stepped::<false>(6);
    let narrow = Wrap { inner: 1u8, len: 2, cap: 8 };
    let broad = Wrap { inner: "w", len: 3, cap: 9 };
    let _ = narrow.slack() + broad.slack();
    let _ = mix("m") + mix([0u8; 5]);
    let _ = vowels("aei");
    let _ = vowels("xyz");
    let _ = largest(&[1u8, 9, 3], 0);
    let _ = largest(&[1.5f32, 0.5], 0.0);
    let _ = padded_len("ab");
    let _ = padded_len([0u8; 2]);
    let _ = spread("spread");
    let _ = spread(vec![7u8]);
    let _ = trailing_spaces("a  ");
    let _ = trailing_spaces("b   ");
    let _ = weigh("w");
    let _ = weigh([1u8]);
    let _ = map_in_part(&mut count, b"map") ^ map_in_part(&mut last, b"map");
    let _ = map_captures_generic(&mut count, b"k", b"map") ^ map_captures_generic(&mut last, &vec![7u8], b"map");
    let _ = map_body_names_param(&mut count, b"map") ^ map_body_names_param(&mut last, b"map");
    let sealed = Sealed { token: 3u8, table: [1, 2, 3, 4, 5, 6, 7, 8], seed: 5, word: 9 };
    let other = Sealed { token: "t", table: [8, 7, 6, 5, 4, 3, 2, 1], seed: 2, word: 4 };
    let coerced: &u32 = &other;
    let _ = sealed.leading_zeros() + coerced;
    let mut rec = Record { name: b"rec".to_vec(), count: 0, total: 0, sent: 0, flags: 0 };
    recount(&mut count, &mut rec, 2);
    recount(&mut last, &mut rec, 3);
    restamp(&mut count, &mut rec, 4);
    restamp(&mut last, &mut rec, 5);
    let _ = tail_after_refused_head(&mut count, &mut rec, 2, b"even");
    let _ = tail_after_refused_head(&mut last, &mut rec, 3, b"odd");
    let mut cur = Cursor { data: vec![1, 2, 3], pos: 0, reads: 0, sum: 0 };
    read_into(&mut count, &mut cur, 2);
    read_into(&mut last, &mut cur, 9);
    {
        let mut v: Vec<u8> = Vec::new();
        let mut k = 0u32;
        read_scanned(&mut v, &mut cur, 2);
        read_scanned(&mut k, &mut cur, 9);
        last.0 ^= v.len() as u8 ^ k as u8;
    }
    {
        let mut a = Cursor { data: vec![1, 2, 3], pos: 0, reads: 0, sum: 0 };
        let mut b = Cursor { data: vec![4, 5], pos: 0, reads: 0, sum: 0 };
        let mut keep: Vec<&[u8]> = Vec::new();
        read_keep(&mut count, &mut a, 2, &mut keep);
        read_keep(&mut last, &mut b, 9, &mut keep);
        last.0 ^= keep.len() as u8;
    }
    {
        let mut a = Cursor { data: vec![1, 2, 3], pos: 0, reads: 0, sum: 0 };
        let mut b = Cursor { data: vec![4, 5], pos: 0, reads: 0, sum: 0 };
        let mut kept = Kept { last: &[] };
        read_store(&mut count, &mut a, 2, &mut kept);
        read_store(&mut last, &mut b, 9, &mut kept);
        last.0 ^= kept.last.len() as u8;
    }
    let boxed = Box::new(Cursor { data: vec![1, 2, 3], pos: 0, reads: 0, sum: 0 });
    let boxed = peek_boxed(&mut count, boxed, 2);
    let _ = peek_boxed(&mut last, boxed, 9);
    let mut src = Src { buf: vec![1, 2, 3], pos: 0 };
    let mut parser = Parser { src: &mut src, depth: 1, sum: 0 };
    next_token(&mut count, &mut parser, 2);
    next_token(&mut last, &mut parser, 9);
    refill(&mut count, &mut cur, 2);
    refill(&mut last, &mut cur, 9);
    refill_twice(&mut count, &mut cur, 2);
    refill_twice(&mut last, &mut cur, 9);
    let _ = refill_total(&mut count, &mut cur, 2) + refill_total(&mut last, &mut cur, 9);
    let lock = std::sync::Mutex::new(5u32);
    let _ = settle(&mut count, &lock, 1) + settle(&mut last, &lock, 2);
    let _ = label(&mut count, 3) + label(&mut last, 4);
    let _ = unlabel(&mut count, 5, true) + unlabel(&mut last, 6, false);
    let _ = handoff(&mut count, 7) + handoff(&mut last, 8);
    own_view(&mut count, 3);
    own_view(&mut last, 5);
    own_iov(&mut count, 3);
    own_iov(&mut last, 5);
    own_non_null(&mut count, 3);
    own_non_null(&mut last, 5);
    own_waker(&mut count, 3);
    own_waker(&mut last, 5);
    own_extent(&mut count, 3);
    own_extent(&mut last, 5);
    parked_view(&mut count, 3);
    parked_view(&mut last, 5);
    rc_parked(&mut count, 3);
    rc_parked(&mut last, 5);
    single_arg_park(&mut count, 3);
    single_arg_park(&mut last, 5);
    tuple_arg_park(&mut count, 3);
    tuple_arg_park(&mut last, 5);
    let _ = trim_stamp(&mut count, 7, 9) + trim_stamp(&mut last, 1, 0);
    let header = FrameHeader { length: 3, kind: 1, flags: 0, stream: 9 };
    let _ = header.write(&mut count) + header.write(&mut last);
    let mut stats = Stats { dots: 0, longest: 0 };
    print_stamp(&mut count, &mut stats, 12, 5);
    print_stamp(&mut last, &mut stats, 3, 0);
    let _ = settings_for((), 3).width + settings_for(7u8, 5).height + settings_for("m", 9).depth;
    let mut printer = Printer { col: 0, out: Vec::new(), errors: 0 };
    let _ = printer.write_str("abc");
    let _ = printer.write_str(b"abc");
    let _ = printer.write_str(vec![1u8, 2]);
    let _small_ledger = Ledger { tag: 1u8, counts: [1, 2, 3, 4, 5, 6], total: 0 };
    let _wide_ledger = Ledger { tag: "wide", counts: [6, 5, 4, 3, 2, 1], total: 0 };
    let _ = macro_checksum("m");
    let _ = macro_checksum([0u8; 3]);
    let mut plain = Conn::<false> { id: 3, open: true };
    let mut secure = Conn::<true> { id: 5, open: true };
    let _ = plain.close() + secure.close();
    let _ = FOLD_TWO + FOLD_THREE;
    let _ = fold_evens(&mut count, b"even") + fold_evens(&mut last, b"odd");
    let rows = [Row { a: 1, b: 2, c: 3 }];
    let _ = hash_row(&mut count, &rows, 0) + hash_row(&mut last, &rows, 1);
    let tagged = Canvas { tag: 1u8, dims: Dims { w: 1, h: 2 } };
    let titled = Canvas { tag: "t", dims: Dims { w: 3, h: 4 } };
    let _ = tagged.area(&mut count) + titled.area(&mut last);
}
