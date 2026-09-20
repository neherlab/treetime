// Several independent bool fields are 2^n representable states.

// Flagged: five bools is 32 states.
struct ShareVerdict {
    share_bp: u64,
    fair: bool,
    rr: bool,
    census: bool,
    stacks: bool,
    counters: bool,
}

// Flagged: three bools is 8 states, and the threshold is three.
struct Refusals {
    too_big: bool,
    exhausted: bool,
    fragmented: bool,
    slots_leaked: u64,
}

// Fine: two bools is under the threshold.
struct Pair {
    a: bool,
    b: bool,
}

// Fine: an explicit repr means the layout is fixed outside Rust, so all eight
// states may genuinely be reachable.
#[repr(C)]
struct HwFlags {
    enabled: bool,
    pending: bool,
    masked: bool,
}

// Fine: a tuple struct has no field names to report.
struct Positional(bool, bool, bool);

fn main() {
    let v = ShareVerdict {
        share_bp: 5000,
        fair: true,
        rr: true,
        census: true,
        stacks: true,
        counters: true,
    };
    println!(
        "{} {} {} {} {} {}",
        v.share_bp, v.fair, v.rr, v.census, v.stacks, v.counters
    );

    let r = Refusals {
        too_big: false,
        exhausted: false,
        fragmented: true,
        slots_leaked: 0,
    };
    println!(
        "{} {} {} {}",
        r.too_big, r.exhausted, r.fragmented, r.slots_leaked
    );

    let p = Pair { a: true, b: false };
    println!("{} {}", p.a, p.b);

    let h = HwFlags {
        enabled: true,
        pending: false,
        masked: false,
    };
    println!("{} {} {}", h.enabled, h.pending, h.masked);

    let t = Positional(true, false, true);
    println!("{} {} {}", t.0, t.1, t.2);
}
