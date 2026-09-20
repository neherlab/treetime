// Names that claim units must agree before they meet in + - or comparison.

struct Timing {
    timeout_ms: u64,
    deadline_ns: u64,
    budget_ms: u64,
    size_bytes: u64,
}

fn mixed_add_is_flagged(t: &Timing) -> u64 {
    t.timeout_ms + t.deadline_ns
}

fn mixed_compare_is_flagged(t: &Timing) -> bool {
    t.timeout_ms < t.deadline_ns
}

fn time_vs_size_is_flagged(t: &Timing) -> bool {
    t.timeout_ms == t.size_bytes
}

fn same_unit_is_fine(t: &Timing) -> u64 {
    t.timeout_ms + t.budget_ms
}

fn alias_units_are_fine(elapsed_millis: u64, budget_ms: u64) -> u64 {
    elapsed_millis + budget_ms
}

fn conversion_is_fine(t: &Timing) -> u64 {
    t.timeout_ms * 1_000_000
}

fn unnamed_operand_is_fine(t: &Timing) -> u64 {
    t.timeout_ms + 5
}

fn method_name_is_flagged(t: &Timing, now_ns: impl Fn() -> u64) -> u64 {
    let _ = &now_ns;
    t.timeout_ms + elapsed_ns()
}

fn elapsed_ns() -> u64 {
    1
}

fn main() {
    let t = Timing {
        timeout_ms: 1,
        deadline_ns: 2,
        budget_ms: 3,
        size_bytes: 4,
    };
    let _ = mixed_add_is_flagged(&t);
    let _ = mixed_compare_is_flagged(&t);
    let _ = time_vs_size_is_flagged(&t);
    let _ = same_unit_is_fine(&t);
    let _ = alias_units_are_fine(5, 6);
    let _ = conversion_is_fine(&t);
    let _ = unnamed_operand_is_fine(&t);
    let _ = method_name_is_flagged(&t, || 7);
}
