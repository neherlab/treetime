// Wildcard arms over small crate-local enums absorb future variants.

enum Op {
    Add,
    Sub,
    Mul,
}

fn wild_is_flagged(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        _ => 0,
    }
}

fn binding_is_flagged(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        other => cost(other),
    }
}

fn cost(_o: Op) -> i32 {
    2
}

fn exhaustive_is_fine(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        Op::Sub => 2,
        Op::Mul => 3,
    }
}

// Extractors ask "is it this one shape?"; future variants are correctly not
// that shape, so `_ => None` and `_ => false` stay silent.
fn extractor_none_is_fine(o: &Op) -> Option<i32> {
    match o {
        Op::Add => Some(1),
        _ => None,
    }
}

fn extractor_false_is_fine(o: &Op) -> bool {
    match o {
        Op::Add => true,
        _ => false,
    }
}

fn extractor_empty_slice_is_fine(o: &Op) -> &'static [u32] {
    match o {
        Op::Add => &[1],
        _ => &[],
    }
}

fn extractor_empty_vec_is_fine(o: &Op) -> Vec<u32> {
    match o {
        Op::Add => vec![1],
        _ => Vec::new(),
    }
}

fn extractor_return_none_is_fine(o: &Op) -> Option<u32> {
    match o {
        Op::Add => Some(1),
        _ => return None,
    }
}

// An #[allow] on the arm itself is honored.
fn arm_allow_is_fine(o: &Op) -> u32 {
    match o {
        Op::Add => 1,
        #[allow(unknown_lints, wildcard_over_own_enum)]
        _ => 0,
    }
}

// `#[non_exhaustive]` constrains downstream crates only. In the defining crate
// — the only place this lint can see a match — rustc still checks
// exhaustiveness, so a wildcard here absorbs future variants exactly as it does
// over an ordinary enum, and on an enum annotated because new variants are
// expected.
#[non_exhaustive]
pub enum Entry {
    Cpu,
    Memory,
    Irq,
}

fn non_exhaustive_wild_is_flagged(e: Entry) -> i32 {
    match e {
        Entry::Cpu => 1,
        _ => 0,
    }
}

// Listing the variants is the fix, and it compiles in the defining crate.
fn non_exhaustive_listed_is_fine(e: Entry) -> i32 {
    match e {
        Entry::Cpu => 1,
        Entry::Memory => 2,
        Entry::Irq => 3,
    }
}

fn foreign_enum_is_fine(x: Option<u32>) -> u32 {
    match x {
        Some(v) => v,
        _ => 0,
    }
}

fn non_enum_is_fine(x: u32) -> u32 {
    match x {
        1 => 2,
        _ => 0,
    }
}

// `b""` and `""` are the empty-slice answer spelled as a literal.
fn extractor_empty_bytestr_is_fine(o: &Op) -> &'static [u8] {
    match o {
        Op::Add => b"+",
        _ => b"",
    }
}

fn extractor_empty_str_is_fine(o: &Op) -> &'static str {
    match o {
        Op::Add => "+",
        _ => "",
    }
}

// The same answers inside braces.
fn extractor_braced_return_is_fine(o: &Op) -> Option<u32> {
    let n = match o {
        Op::Add => 1,
        _ => {
            return None;
        }
    };
    Some(n)
}

fn extractor_braced_tail_is_fine(o: &Op) -> Option<u32> {
    match o {
        Op::Add => Some(1),
        _ => {
            // Nothing else has a cost.
            None
        }
    }
}

// A statement before the answer is behavior the arm gives future variants.
fn braced_statement_is_flagged(o: &Op) -> Option<u32> {
    match o {
        Op::Add => Some(1),
        _ => {
            note();
            None
        }
    }
}

fn note() {}

// A binding arm whose whole body is another match on the binding hands the
// dispatch on; the inner match is judged on its own.
fn redispatch_is_fine(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        other => match other {
            Op::Add => 0,
            Op::Sub => 2,
            Op::Mul => 3,
        },
    }
}

macro_rules! dispatch {
    ($o:expr) => {
        match $o {
            Op::Add => 1,
            Op::Sub => 2,
            Op::Mul => 3,
        }
    };
}

fn redispatch_via_macro_is_fine(o: Op) -> i32 {
    match o {
        Op::Mul => 30,
        other => dispatch!(other),
    }
}

// Only the inner match's own wildcard is flagged.
fn redispatch_flags_inner_only(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        other => match other {
            Op::Sub => 2,
            _ => 0,
        },
    }
}

// A statement before the inner match runs for future variants too.
fn redispatch_after_statement_is_flagged(o: Op) -> i32 {
    match o {
        Op::Add => 1,
        other => {
            note();
            match other {
                Op::Add => 0,
                Op::Sub => 2,
                Op::Mul => 3,
            }
        }
    }
}

fn main() {
    let _ = wild_is_flagged(Op::Add);
    let _ = binding_is_flagged(Op::Sub);
    let _ = exhaustive_is_fine(Op::Mul);
    let _ = extractor_none_is_fine(&Op::Add);
    let _ = extractor_false_is_fine(&Op::Sub);
    let _ = extractor_empty_slice_is_fine(&Op::Add);
    let _ = extractor_empty_vec_is_fine(&Op::Sub);
    let _ = extractor_return_none_is_fine(&Op::Mul);
    let _ = arm_allow_is_fine(&Op::Add);
    let _ = non_exhaustive_wild_is_flagged(Entry::Cpu);
    let _ = non_exhaustive_listed_is_fine(Entry::Memory);
    let _ = non_exhaustive_listed_is_fine(Entry::Irq);
    let _ = foreign_enum_is_fine(Some(1));
    let _ = non_enum_is_fine(1);
    let _ = extractor_empty_bytestr_is_fine(&Op::Add);
    let _ = extractor_empty_str_is_fine(&Op::Sub);
    let _ = extractor_braced_return_is_fine(&Op::Mul);
    let _ = extractor_braced_tail_is_fine(&Op::Add);
    let _ = braced_statement_is_flagged(&Op::Sub);
    let _ = redispatch_is_fine(Op::Mul);
    let _ = redispatch_via_macro_is_fine(Op::Add);
    let _ = redispatch_flags_inner_only(Op::Sub);
    let _ = redispatch_after_statement_is_flagged(Op::Mul);
}
