// FAILING: suppressions without a reason.
#[allow(dead_code)] // should warn: unjustified_suppression
struct Unused;

fn body_allow() {
    #[allow(unused_variables)] // should warn: unjustified_suppression
    let x = 1;
}

// PASSING: suppressions that carry a reason (and prefer `#[expect]`).
#[expect(dead_code, reason = "kept for the public API surface")]
struct KeptForApi;

#[allow(dead_code, reason = "constructed only through FFI")]
struct FromFfi;

fn main() {}
