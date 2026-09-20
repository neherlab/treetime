// Panics rustc lowers to a MIR `Assert` terminator, one root per family. None
// of these is a call in HIR at any spelling, so before the MIR walk every one
// of them was invisible to a `never` list; `forbidden_reach.rs`'s `index_root`
// covers the bounds-check family the same way.

// Flagged: `a + b` on integers becomes `Assert(Overflow(Add, ..))`, whose
// panic is `core::panicking::panic_const::panic_const_add_overflow`. This is
// the one family here that depends on the profile: with `-C overflow-checks`
// off there is no assert to see, and the ui harness runs with them on.
fn add_overflow_root(a: u32, b: u32) -> u32 {
    a + b
}

// Flagged: `Assert(DivisionByZero(..))`, emitted whatever the profile asks
// for, unlike the overflow family above.
fn div_zero_root(a: u32, b: u32) -> u32 {
    a / b
}

// Flagged: `Assert(RemainderByZero(..))`, the sibling of the division check
// and a separate lang item, so a `never` list naming one does not catch the
// other.
fn rem_zero_root(a: u32, b: u32) -> u32 {
    a % b
}

// Flagged: unary negation of a signed integer is `Assert(OverflowNeg(..))`,
// which reaches a third overflow lang item again distinct from `Add`'s.
fn neg_overflow_root(a: i32) -> i32 {
    -a
}

// Fine: a root banned from the add-overflow panic whose body asserts twice,
// on division and remainder by zero. The control for the four above -- it
// proves a finding tracks which assert the body has and not merely that it
// has one. Comparing the two results rather than adding them is the whole
// point: `a / b + a % b` reaches `panic_const_add_overflow` through the `+`,
// so it is a fifth positive and not a control at all.
fn wrong_family_root(a: u32, b: u32) -> u32 {
    if a / b > a % b { a } else { b }
}

// Fine: `wrapping_add` compiles to no assert at all, so a declared root doing
// only that reaches nothing. The control for "the ban is live but silent".
fn no_assert_root(a: u32, b: u32) -> u32 {
    a.wrapping_add(b)
}

fn main() {
    let _ = add_overflow_root(1, 2);
    let _ = div_zero_root(6, 3);
    let _ = rem_zero_root(6, 3);
    let _ = neg_overflow_root(1);
    let _ = wrong_family_root(6, 3);
    let _ = no_assert_root(1, 2);
}
