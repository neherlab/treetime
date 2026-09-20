// Each family of lints is a lint group named `mordant_<family>`; allowing
// the group silences its members.

enum Op {
    Add,
    Sub,
    Mul,
}

// Flagged: `wildcard_over_own_enum`, nothing allows it.
fn wild(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        _ => 0,
    }
}

// Silenced: `wildcard_over_own_enum` is in `enums`.
#[allow(mordant_enums)]
fn wild_allowed(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        _ => 0,
    }
}

// Still flagged: `discarded_error` is in `errors`, not `enums`.
#[allow(mordant_enums)]
fn other_family() {
    "1".parse::<u8>().ok();
}

fn main() {}
