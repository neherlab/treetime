// compile-flags: --edition=2024
//
// A `Some` (or `Ok`) taken only under a further condition, with the value
// that fails it dropped into the same handling as `None`: the condition
// belongs where the value is produced, not at this consumer.

struct Job {
    ready: bool,
    id: u32,
}

fn next() -> Option<Job> {
    Some(Job { ready: true, id: 1 })
}

fn parse() -> Result<u32, String> {
    Ok(3)
}

fn wait() -> u32 {
    0
}

// Flagged: the guard reads `j`, and a `Some` failing it falls to `_` with
// `None`.
fn guarded_arm_with_wildcard() -> u32 {
    match next() {
        Some(j) if j.ready => j.id,
        _ => wait(),
    }
}

// Flagged: the same fallthrough spelled as two arms with identical bodies,
// whatever their order.
fn guarded_arm_with_split_rest() -> u32 {
    match next() {
        None => wait(),
        Some(j) if j.ready => j.id,
        Some(_) => wait(),
    }
}

// Flagged: `Some(_) | None` is `_` written out.
fn guarded_arm_with_or_rest() -> u32 {
    match next() {
        Some(j) if j.ready && j.id > 0 => j.id,
        Some(_) | None => wait(),
    }
}

// Flagged: an `Ok` failing the guard is handled as the error is.
fn guarded_ok_arm() -> u32 {
    match parse() {
        Ok(n) if n > 2 => n,
        Err(_) | Ok(_) => wait(),
    }
}

// Flagged: the let-chain spelling.
fn let_chain_with_else() -> u32 {
    if let Some(j) = next()
        && j.ready
    {
        j.id
    } else {
        wait()
    }
}

// Flagged: the condition as an inner `if` whose `else` repeats the outer.
fn nested_if_repeating_else() -> u32 {
    if let Some(j) = next() {
        if j.ready { j.id } else { wait() + 1 }
    } else {
        wait() + 1
    }
}

// Fine: the guard never reads the binding; that is a flag beside the
// option, a different smell.
fn guard_on_other_state(open: bool) -> u32 {
    match next() {
        Some(j) if open => j.id,
        _ => wait(),
    }
}

// Fine: the failing `Some` and `None` are handled by different code.
fn rest_arms_differ() -> u32 {
    match next() {
        Some(j) if j.ready => j.id,
        Some(_) => wait() + 1,
        None => wait(),
    }
}

// Fine: the inner `else` and the outer one differ.
fn nested_if_different_else() -> u32 {
    if let Some(j) = next() {
        if j.ready { j.id } else { wait() + 1 }
    } else {
        wait()
    }
}

// Fine: a condition, not handling; this is what the fixed consumer reads.
fn is_some_and_condition() -> u32 {
    if next().is_some_and(|j| j.ready) {
        1
    } else {
        wait()
    }
}

// Fine: `matches!` with a guard is a condition too.
fn matches_condition() -> u32 {
    if matches!(next(), Some(j) if j.ready) {
        1
    } else {
        wait()
    }
}

enum Slot {
    Full(Job),
    Empty,
}

// Fine: not an `Option` or a `Result`.
fn local_enum_with_guard(s: Slot) -> u32 {
    match s {
        Slot::Full(j) if j.ready => j.id,
        Slot::Full(_) | Slot::Empty => wait(),
    }
}

// Fine: the second `Some` arm does its own work with the value.
fn second_some_arm_works() -> u32 {
    match next() {
        Some(j) if j.ready => j.id,
        Some(j) => j.id + 100,
        None => wait(),
    }
}

macro_rules! ready_or_wait {
    ($e:expr) => {
        match $e {
            Some(j) if j.ready => j.id,
            _ => wait(),
        }
    };
}

// Fine: inside a macro expansion.
fn inside_macro() -> u32 {
    ready_or_wait!(next())
}

// Fine: no `else`, so nothing handles either case here.
fn let_chain_without_else() -> u32 {
    if let Some(j) = next()
        && j.ready
    {
        return j.id;
    }
    2
}

// Fine: the fallback does nothing.
fn empty_fallback(total: &mut u32) {
    match next() {
        Some(j) if j.ready => *total += j.id,
        _ => {}
    }
}

fn main() {
    let _ = guarded_arm_with_wildcard();
    let _ = guarded_arm_with_split_rest();
    let _ = guarded_arm_with_or_rest();
    let _ = guarded_ok_arm();
    let _ = let_chain_with_else();
    let _ = nested_if_repeating_else();
    let _ = guard_on_other_state(true);
    let _ = rest_arms_differ();
    let _ = nested_if_different_else();
    let _ = is_some_and_condition();
    let _ = matches_condition();
    let _ = local_enum_with_guard(Slot::Empty);
    let _ = local_enum_with_guard(Slot::Full(Job { ready: false, id: 0 }));
    let _ = second_some_arm_works();
    let _ = inside_macro();
    let _ = let_chain_without_else();
    empty_fallback(&mut 0);
}
