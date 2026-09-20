// The same `match` over one enum written out arm for arm in two functions is
// a method the enum lacks; matches that differ anywhere are different code.

use std::fmt;

pub enum Step {
    Open,
    Read,
    Parse,
}

pub struct Failure {
    pub step: Step,
    pub code: u32,
}

pub fn describe(f: &Failure) -> (u32, &'static str) {
    match f.step {
        Step::Open => (f.code, "open"),
        Step::Read => (f.code + 1, "read"),
        Step::Parse => (f.code + 2, "parse"),
    }
}

// Flagged: the same three arms as `describe`, reading a local of the same
// type under another name.
pub fn describe_again(cause: &Failure, loud: bool) -> Option<(u32, &'static str)> {
    if !loud {
        return None;
    }
    Some(match cause.step {
        Step::Open => (cause.code, "open"),
        Step::Read => (cause.code + 1, "read"),
        Step::Parse => (cause.code + 2, "parse"),
    })
}

// Flagged once, against `describe`: a third copy does not also pair with the
// second.
pub fn describe_third(f: &Failure) -> (u32, &'static str) {
    let d = match f.step {
        Step::Open => (f.code, "open"),
        Step::Read => (f.code + 1, "read"),
        Step::Parse => (f.code + 2, "parse"),
    };
    d
}

// Fine: one arm body differs.
pub fn describe_short(f: &Failure) -> (u32, &'static str) {
    match f.step {
        Step::Open => (f.code, "open"),
        Step::Read => (f.code + 1, "read"),
        Step::Parse => (0, "parse"),
    }
}

pub fn report(f: &Failure) -> String {
    match f.step {
        Step::Open => format!("failed to open {}", f.code),
        Step::Read => format!("failed to read {}", f.code),
        Step::Parse => String::from("failed to parse"),
    }
}

// Fine, though a copy: arguments to a macro compare by their tokens, so the
// renamed local inside `format!` makes these arms different text.
pub fn report_renamed(cause: &Failure) -> String {
    match cause.step {
        Step::Open => format!("failed to open {}", cause.code),
        Step::Read => format!("failed to read {}", cause.code),
        Step::Parse => String::from("failed to parse"),
    }
}

#[derive(Clone, Copy)]
pub enum Level {
    Low,
    Mid,
    High,
}

pub struct Gauge {
    pub level: Level,
    pub scale: u32,
}

pub struct Meter {
    pub level: Level,
    pub scale: u32,
    pub live: bool,
}

pub fn gauge_units(g: &Gauge) -> u64 {
    match g.level {
        Level::Low => 1,
        Level::Mid => u64::from(g.scale),
        Level::High => 100,
    }
}

// Fine: the arms read `m`, whose type is not `g`'s, so this is not the same
// code even though it is spelled alike.
pub fn meter_units(m: &Meter) -> u64 {
    match m.level {
        Level::Low => 1,
        Level::Mid => u64::from(m.scale),
        Level::High => 100,
    }
}

pub fn low_offset(l: Level) -> Option<u32> {
    match l {
        Level::Low => Some(3),
        _ => None,
    }
}

// Fine: one counted arm and a catch-all is a test, not a table.
pub fn low_offset_again(l: Level) -> Option<u32> {
    match l {
        Level::Low => Some(3),
        _ => None,
    }
}

// Fine: `Display` and `Debug` are the enum's own trait impls, which have to
// match every variant to exist.
impl fmt::Display for Level {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Level::Low => f.write_str("low"),
            Level::Mid => f.write_str("mid"),
            Level::High => f.write_str("high"),
        }
    }
}

impl fmt::Debug for Level {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("Level::")?;
        match self {
            Level::Low => f.write_str("low"),
            Level::Mid => f.write_str("mid"),
            Level::High => f.write_str("high"),
        }
    }
}

pub fn some_or_zero(o: Option<u32>) -> u32 {
    match o {
        Some(v) => v + 1,
        None => 0,
    }
}

// Fine: `Option` is a standard-library enum; a repeated match on it is an
// idiom, not a table the crate could turn into a method.
pub fn some_or_zero_again(o: Option<u32>) -> u32 {
    match o {
        Some(v) => v + 1,
        None => 0,
    }
}

pub enum Shape {
    Dot,
    Line(u32),
    Rect(u32, u32),
}

pub fn area(s: &Shape) -> u32 {
    match s {
        Shape::Dot => 0,
        Shape::Line(_) => 0,
        Shape::Rect(w, h) => match (w, h) {
            (0, h) => *h,
            (w, h) => w * h,
        },
    }
}

// Flagged once: the outer match repeats `area`'s; the tuple match inside it
// is covered by that report and not named again.
pub fn area_again(shape: &Shape) -> u32 {
    let cells = match shape {
        Shape::Dot => 0,
        Shape::Line(_) => 0,
        Shape::Rect(w, h) => match (w, h) {
            (0, h) => *h,
            (w, h) => w * h,
        },
    };
    cells.next_power_of_two()
}

// Fine: reads two free locals where `area` reads one.
pub fn area_scaled(s: &Shape, k: u32) -> u32 {
    match s {
        Shape::Dot => 0,
        Shape::Line(_) => k,
        Shape::Rect(w, h) => w * h,
    }
}

#[derive(Clone, Copy)]
pub enum Axis {
    Across,
    Down,
}

// Flagged: the same match a second time in one body, over the same locals.
pub fn reach(axis: Axis, w: u32, h: u32) -> (u32, u32) {
    let near = match axis {
        Axis::Across => w,
        Axis::Down => h + 1,
    };
    let far = match axis {
        Axis::Across => w,
        Axis::Down => h + 1,
    };
    (near, far * 2)
}

// Fine: two matches in one body reading the same locals in different
// places. `w` stands where `h` does in the first arm and where `w` itself
// does in the second, so no renaming of locals turns one into the other.
pub fn corners(axis: Axis, w: u32, h: u32, scale: u32) -> ((u32, u32), (u32, u32)) {
    let p = match axis {
        Axis::Across => (w, scale),
        Axis::Down => (scale, w),
    };
    let q = match axis {
        Axis::Across => (h, w),
        Axis::Down => (w, w),
    };
    (p, q)
}

// Fine: a macro that expands its argument twice makes two matches out of
// one piece of source, which is still one copy.
macro_rules! both_ways {
    ($e:expr) => {
        ($e, $e)
    };
}

pub fn level_pair(l: Level) -> (u8, u8) {
    both_ways!(match l {
        Level::Low => 0,
        Level::Mid => 5,
        Level::High => 9,
    })
}

fn main() {}
