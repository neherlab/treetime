#![allow(dead_code, unused_variables)]

pub const MAX_WORDS: u32 = 1024;
pub const MAX_TYPES: u32 = 256;
pub const MAX_BLOCKS: u32 = 64;

pub enum Ceiling {
    Words,
    Types,
    Blocks,
}

/// The shape: `limit` restates the constant `ceiling` already names, and the
/// pairing is held by three call sites rather than by a type.
pub struct Exceeded {
    ceiling: Ceiling,
    limit: u32,
    wanted: u32,
}

pub fn words(wanted: u32) -> Exceeded {
    Exceeded {
        ceiling: Ceiling::Words,
        limit: MAX_WORDS,
        wanted,
    }
}

pub fn types(wanted: u32) -> Exceeded {
    Exceeded {
        ceiling: Ceiling::Types,
        limit: MAX_TYPES,
        wanted,
    }
}

pub fn blocks(wanted: u32) -> Exceeded {
    Exceeded {
        ceiling: Ceiling::Blocks,
        limit: MAX_BLOCKS,
        wanted,
    }
}

// ---------------------------------------------------------------------------
// The `Option`-beside-an-enum shape: `remaining` is meaningful for exactly one
// variant of `header`, and the two are written independently.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
pub enum Header {
    Narrow,
    Wide,
}

pub struct Entries<'a> {
    buf: &'a [u8],
    header: Header,
    remaining: Option<u32>,
}

pub fn narrow(buf: &[u8]) -> Entries<'_> {
    Entries {
        buf,
        header: Header::Narrow,
        remaining: None,
    }
}

pub fn wide(buf: &[u8], count: u32) -> Entries<'_> {
    Entries {
        buf,
        header: Header::Wide,
        remaining: Some(count),
    }
}

// ---------------------------------------------------------------------------
// Must not fire.
// ---------------------------------------------------------------------------

/// Two columns of bare literals. Any two test-table rows that differ are in
/// bijection; neither field names a value a reader could have got wrong.
pub struct Row {
    tsc: u64,
    duration: u64,
}

pub fn row_a() -> Row {
    Row {
        tsc: 0x0123,
        duration: 4,
    }
}

pub fn row_b() -> Row {
    Row {
        tsc: 0x4567,
        duration: 5,
    }
}

/// One construction site cannot exhibit a correspondence. `b` and `c` are
/// derived from `a` on purpose so `intact` can detect corruption — the
/// redundancy is the oracle, and a lint that fired here would be telling a
/// test to delete its own check.
#[derive(PartialEq)]
pub struct Tri {
    a: u64,
    b: u64,
    c: u64,
}

impl Tri {
    pub fn of(a: u64) -> Self {
        Self {
            a,
            b: a.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            c: a ^ 0xDEAD_BEEF_CAFE_F00D,
        }
    }
    pub fn intact(&self) -> bool {
        *self == Self::of(self.a)
    }
}

/// A pairing that is not one-for-one: `Ceiling::Words` appears beside two
/// different limits, so `limit` is not decided by `ceiling`.
pub struct Loose {
    ceiling: Ceiling,
    limit: u32,
}

pub fn loose_a() -> Loose {
    Loose {
        ceiling: Ceiling::Words,
        limit: MAX_WORDS,
    }
}

pub fn loose_b() -> Loose {
    Loose {
        ceiling: Ceiling::Words,
        limit: MAX_TYPES,
    }
}

pub fn loose_c() -> Loose {
    Loose {
        ceiling: Ceiling::Types,
        limit: MAX_BLOCKS,
    }
}

/// A two-row table. Every column is in bijection with every other, because
/// two rows that differ anywhere differ everywhere; `budget` names a constant
/// but the constant is this row's own datum, not a restatement of `name`.
/// Measured firing four times on one such table before `decides` narrowed.
pub struct Entry {
    name: &'static str,
    why: &'static str,
    budget: u32,
}

pub const FAST: u32 = 20;
pub const SLOW: u32 = 90;

pub fn table() -> [Entry; 2] {
    [
        Entry {
            name: "ctx_switch",
            why: "the second headline number",
            budget: FAST,
        },
        Entry {
            name: "wake_latency",
            why: "the first headline number",
            budget: SLOW,
        },
    ]
}

/// An explicit `repr` means the layout is dictated from outside Rust, so a
/// field restating what a sibling implies may be the wire format's doing.
#[repr(C)]
pub struct Wire {
    kind: u32,
    width: u32,
}

pub fn wire_a() -> Wire {
    Wire { kind: 1, width: 4 }
}

pub fn wire_b() -> Wire {
    Wire { kind: 2, width: 8 }
}

/// A `..base` literal overrides part of a value that already exists, so the
/// fields written are the ones deliberately being made to differ.
#[derive(Clone, Copy)]
pub struct Record {
    turn: u32,
    steps: u32,
}

pub const BLANK: Record = Record { turn: 0, steps: 0 };

pub fn records() -> [Record; 2] {
    [
        Record { turn: 1, ..BLANK },
        Record {
            turn: 2,
            steps: 9,
            ..BLANK
        },
    ]
}

/// `owner` is filled in after construction, so the literals are not the only
/// thing deciding it: every `Detached` site starts it at `None`, and `attach`
/// pairs it with `Detached` anyway.
#[derive(Clone, Copy, PartialEq)]
pub enum Source {
    Attached,
    Detached,
}

pub struct Job {
    source: Source,
    owner: Option<u32>,
}

pub fn attached(owner: u32) -> Job {
    Job {
        source: Source::Attached,
        owner: Some(owner),
    }
}

pub fn detached() -> Job {
    Job {
        source: Source::Detached,
        owner: None,
    }
}

pub fn attach(job: &mut Job, owner: u32) {
    job.owner = Some(owner);
}

fn main() {}
