// A private enum variant that is constructed but never named by a pattern
// outside the enum's own impls carries structure nobody reads.

use std::fmt;

#[derive(Debug)]
enum LoadError {
    NotFound,
    Corrupt(String),
}

impl fmt::Display for LoadError {
    // The enum's own impls must match every variant; this does not count as
    // the crate reading the structure.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LoadError::NotFound => write!(f, "not found"),
            LoadError::Corrupt(why) => write!(f, "corrupt: {why}"),
        }
    }
}

fn load(x: u32) -> Result<u32, LoadError> {
    match x {
        0 => Err(LoadError::NotFound),
        1 => Err(LoadError::Corrupt("bad header".to_owned())),
        n => Ok(n),
    }
}

fn handle(x: u32) -> u32 {
    match load(x) {
        Ok(n) => n,
        // NotFound is distinguished; Corrupt only ever reaches the catch-all,
        // so its payload is never read and the variant is flagged.
        Err(LoadError::NotFound) => 0,
        Err(other) => {
            eprintln!("{other}");
            0
        }
    }
}

// Every variant is distinguished somewhere; nothing to flag.
enum Poll {
    Ready(u32),
    Pending,
}

fn drain(p: Poll) -> u32 {
    match p {
        Poll::Ready(n) => n,
        Poll::Pending => 0,
    }
}

// Fine: an inherent accessor is the crate genuinely reading the structure,
// unlike a trait impl, so both variants count as named.
enum Ceiling {
    Budget,
    Tripwire(u64),
}

impl Ceiling {
    fn tenths(&self, budget: u64) -> u64 {
        match self {
            Ceiling::Budget => budget,
            Ceiling::Tripwire(t) => *t,
        }
    }
}

fn resolve(budget: u64) -> u64 {
    Ceiling::Budget.tenths(budget) + Ceiling::Tripwire(3).tenths(budget)
}

// `==` against a variant distinguishes it like a pattern: `Off` is fine.
// `Slow` and `Crawl` are constructed but only ever reach the fall-through
// together, so nothing tells them apart and both are flagged.
#[derive(PartialEq)]
enum Mode {
    Fast,
    Slow,
    Crawl,
    Off,
}

fn speed(m: Mode) -> u32 {
    if m == Mode::Off {
        return 0;
    }
    if matches!(m, Mode::Fast) { 2 } else { 1 }
}

fn speeds() -> u32 {
    speed(Mode::Fast) + speed(Mode::Slow) + speed(Mode::Crawl) + speed(Mode::Off)
}

// Fine: with `Active` and `Inactive` both named, `Done` is the only variant
// left, so `!= Inactive` after `!= Active` reaches exactly it.
#[derive(PartialEq)]
enum Status {
    Active,
    Inactive,
    Done,
}

fn step(s: &mut Status) -> bool {
    if *s == Status::Active {
        *s = Status::Done;
        return true;
    }
    if *s != Status::Inactive {
        return false;
    }
    *s = Status::Active;
    true
}

fn run() -> bool {
    let mut s = Status::Inactive;
    step(&mut s) && step(&mut s) && !step(&mut s)
}

// A cast reads every discriminant, so `High` is read even though only `Low`
// is ever named by a pattern.
enum Level {
    Low = 1,
    High = 2,
}

fn code(l: Level) -> u8 {
    l as u8
}

fn is_low(l: &Level) -> bool {
    matches!(l, Level::Low)
}

fn levels() -> u8 {
    u8::from(is_low(&Level::Low)) + code(Level::Low) + code(Level::High)
}

// No pattern anywhere names any variant, so matching is not how this enum is
// consumed and the lint stays silent about all of it.
enum Trace {
    Enter,
    Exit,
}

impl fmt::Display for Trace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Trace::Enter => write!(f, "enter"),
            Trace::Exit => write!(f, "exit"),
        }
    }
}

fn trace() {
    println!("{} {}", Trace::Enter, Trace::Exit);
}

// A tuple constructor handed to `map_err` constructs the variant just as a
// written-out call does: `Malformed` is built that way and never named, so it
// is flagged; `Missing` is built the same way and named, so it is not.
enum ParseError {
    Missing,
    Malformed(std::num::ParseIntError),
}

fn parse(s: &str) -> Result<u32, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Missing);
    }
    s.parse::<u32>().map_err(ParseError::Malformed)
}

fn parse_or_zero(s: &str) -> u32 {
    match parse(s) {
        Ok(n) => n,
        Err(ParseError::Missing) => 0,
        Err(_) => 1,
    }
}

fn main() {
    let _ = handle(0);
    let _ = resolve(10);
    let _ = drain(Poll::Ready(1));
    let _ = drain(Poll::Pending);
    trace();
    let _ = speeds();
    let _ = run();
    let _ = levels();
    let _ = parse_or_zero("7");
}
