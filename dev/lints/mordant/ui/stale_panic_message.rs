// Panic and expect messages must name things that still exist.

struct Table {
    slots: Vec<u32>,
}

fn crate_name_is_fine(t: &Table, i: usize) -> u32 {
    *t.slots.get(i).expect("index checked against `slots` length")
}

fn stale_expect_is_flagged(t: &Table, i: usize) -> u32 {
    *t.slots.get(i).expect("guarded upstream by `frame_lock`")
}

fn stale_macro_is_flagged(x: u32) -> u32 {
    match x {
        0 => 0,
        _ => unreachable!("`normalizer` rejects nonzero values before this"),
    }
}

fn live_macro_is_fine(x: u32) -> u32 {
    if x > 10 {
        panic!("`clamp` should have bounded this");
    }
    clamp(x)
}

fn clamp(x: u32) -> u32 {
    x.min(10)
}

fn plain_message_is_fine(x: u32) -> u32 {
    assert!(x < 100, "value out of range");
    x
}

fn main() {
    let t = Table { slots: vec![1] };
    let _ = crate_name_is_fine(&t, 0);
    let _ = stale_expect_is_flagged(&t, 0);
    let _ = stale_macro_is_flagged(0);
    let _ = live_macro_is_fine(3);
    let _ = plain_message_is_fine(5);
}
