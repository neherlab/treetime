// A string field or local only ever holding one of a closed set of literals, then compared against them.

// Flagged: every store of `kind` is one of two literals and readers compare it.
struct Task {
    kind: &'static str,
    id: u32,
}

fn download(id: u32) -> Task {
    Task { kind: "download", id }
}

fn extract(id: u32) -> Task {
    let mut t = Task { kind: "download", id };
    t.kind = "extract";
    t
}

fn run(t: &Task) -> u32 {
    if t.kind == "download" { t.id } else { 0 }
}

// Flagged: owned bytes built from literals, matched on with literal arms.
struct Phase {
    name: Vec<u8>,
}

fn phases() -> [Phase; 2] {
    [
        Phase {
            name: b"resolve".to_vec(),
        },
        Phase {
            name: b"link".to_vec(),
        },
    ]
}

fn is_link(p: &Phase) -> bool {
    match p.name.as_slice() {
        b"link" => true,
        _ => false,
    }
}

// Flagged: a `String` compared through a comparing method.
struct Mode {
    label: String,
}

fn modes(fast: bool) -> Mode {
    if fast {
        Mode {
            label: String::from("fast"),
        }
    } else {
        Mode {
            label: "slow".to_string(),
        }
    }
}

fn is_fast(m: &Mode) -> bool {
    m.label.eq_ignore_ascii_case("FAST")
}

// Flagged: the struct is public but the field is not, so every store is still in this crate.
pub struct Conn {
    scheme: &'static [u8],
    pub port: u16,
}

pub fn plain() -> Conn {
    Conn {
        scheme: b"ws",
        port: 80,
    }
}

pub fn secure() -> Conn {
    Conn {
        scheme: b"wss",
        port: 443,
    }
}

pub fn default_port(c: &Conn) -> bool {
    matches!(c.scheme, b"wss") == (c.port == 443)
}

// Flagged: a local assigned literals on every path and then only compared.
fn severity(major: bool, minor: bool) -> u8 {
    let mut color = "green";
    if major {
        color = "red";
    } else if minor {
        color = "yellow";
    }
    if color == "red" {
        2
    } else if color == "yellow" {
        1
    } else {
        0
    }
}

// Flagged: one store, but it chooses between two literals.
fn tier(n: u32) -> bool {
    let level = if n > 9 { "hi" } else { "lo" };
    matches!(level, "hi")
}

// Flagged: a derived `Clone` copies the field, which adds no new value.
#[derive(Clone)]
struct Stage {
    step: &'static str,
}

fn stages() -> [Stage; 3] {
    let first = Stage { step: "parse" };
    [first.clone(), Stage { step: "emit" }, first]
}

fn is_parse(s: &Stage) -> bool {
    s.step == "parse"
}

// Fine: the local is only formatted.
fn banner(wide: bool) -> String {
    let rule = if wide { "====" } else { "--" };
    format!("{rule} title {rule}")
}

// Fine: one store comes from a call, so the set is open.
fn from_env(f: fn() -> &'static str) -> bool {
    let mut mode = "auto";
    if f().is_empty() {
        mode = f();
    }
    mode == "auto"
}

// Fine: a `ref mut` binding of the local is a write the lint cannot read.
fn rebound(flip: bool, f: fn() -> &'static str) -> bool {
    let mut mode = "x";
    if flip {
        mode = "y";
    }
    {
        let ref mut r = mode;
        *r = f();
    }
    mode == "x"
}

// Fine: mutated in place through a method taking `&mut self`.
fn grown() -> bool {
    let mut s = String::from("a");
    s.push('b');
    s = "c".to_string();
    s == "ab"
}

// Fine: one store is not a literal, so the set is open.
struct Named {
    name: &'static str,
}

fn named(n: &'static str) -> [Named; 3] {
    [Named { name: "a" }, Named { name: "b" }, Named { name: n }]
}

fn is_a(n: &Named) -> bool {
    n.name == "a"
}

// Fine: only ever formatted, never compared.
struct Label {
    text: &'static str,
}

fn labels() -> [Label; 2] {
    [Label { text: "one" }, Label { text: "two" }]
}

fn show(l: &Label) -> String {
    format!("{}", l.text)
}

// Fine: a single literal is a constant, not a state.
struct Fixed {
    tag: &'static str,
}

fn fixed() -> [Fixed; 2] {
    [Fixed { tag: "x" }, Fixed { tag: "x" }]
}

fn is_x(f: &Fixed) -> bool {
    f.tag == "x"
}

// Fine: the field is mutated in place, so the literals are only seeds.
struct Buf {
    text: String,
}

fn bufs() -> [Buf; 2] {
    [
        Buf {
            text: "a".to_owned(),
        },
        Buf {
            text: "b".to_owned(),
        },
    ]
}

fn grow(b: &mut Buf) -> bool {
    b.text.push('!');
    b.text == "a!"
}

// Fine: part of the value is overwritten in place through an index.
struct Bytes {
    buf: Box<[u8]>,
}

fn bytes() -> [Bytes; 2] {
    [
        Bytes {
            buf: b"ab".to_vec().into_boxed_slice(),
        },
        Bytes {
            buf: b"cd".to_vec().into_boxed_slice(),
        },
    ]
}

fn poke(b: &mut Bytes) -> bool {
    b.buf[0] = b'z';
    &b.buf[..] == b"zb"
}

// Fine: destructured through `&mut`, so the binding is a `&mut` to the field.
struct Job {
    state: String,
    tries: u32,
}

fn jobs() -> [Job; 2] {
    [
        Job {
            state: "queued".into(),
            tries: 0,
        },
        Job {
            state: "done".into(),
            tries: 0,
        },
    ]
}

fn advance(js: &mut [Job], next: &str) -> bool {
    for Job { state, tries } in js.iter_mut() {
        *tries += 1;
        state.clear();
        state.push_str(next);
    }
    js[0].state == "done"
}

// Fine: an explicit `ref mut` in a match arm writes the field.
struct Step {
    kind: &'static str,
}

impl Step {
    fn set(&mut self, s: &'static str) {
        match *self {
            Step { ref mut kind } => *kind = s,
        }
    }
}

fn steps() -> [Step; 2] {
    [Step { kind: "download" }, Step { kind: "extract" }]
}

fn is_download(s: &Step) -> bool {
    s.kind == "download"
}

// Fine: `..base` fills the field from somewhere this site does not spell.
struct Opts {
    level: &'static str,
    n: u32,
}

fn opts(base: &Opts) -> [Opts; 3] {
    [
        Opts { level: "hi", n: 0 },
        Opts { level: "lo", n: 1 },
        Opts { n: 2, ..*base },
    ]
}

fn is_hi(o: &Opts) -> bool {
    o.level == "hi" && o.n < 5
}

// Fine: exported, so other crates may store anything.
pub struct Public {
    pub state: &'static str,
}

pub fn publics() -> [Public; 2] {
    [Public { state: "on" }, Public { state: "off" }]
}

pub fn is_on(p: &Public) -> bool {
    p.state == "on"
}

// Fine: an explicit repr means the layout is fixed from outside.
#[repr(C)]
struct Wire {
    op: &'static str,
}

fn wires() -> [Wire; 2] {
    [Wire { op: "get" }, Wire { op: "set" }]
}

fn is_get(w: &Wire) -> bool {
    w.op == "get"
}

fn main() {
    let _ = run(&download(1)) + run(&extract(2));
    let _ = phases().iter().any(is_link);
    let _ = is_fast(&modes(true));
    let _ = severity(true, false) + tier(3) as u8;
    let _ = stages().iter().any(is_parse);
    let _ = banner(true);
    let _ = from_env(|| "x") || rebound(true, || "z") || grown();
    let _ = named("c").iter().any(is_a);
    let _ = labels().iter().map(show).count();
    let _ = fixed().iter().any(is_x);
    let _ = bufs().iter_mut().any(grow);
    let _ = bytes().iter_mut().any(poke);
    let _ = advance(&mut jobs(), "next");
    let mut all = steps();
    all[0].set("verify");
    let _ = all.iter().any(is_download);
    let base = Opts { level: "mid", n: 9 };
    let _ = opts(&base).iter().any(is_hi);
    let _ = publics().iter().any(is_on);
    let _ = wires().iter().any(is_get);
}
