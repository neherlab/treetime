#![allow(dead_code, unused_variables, clippy::all, reason = "ui fixture exercises calls that are never executed")]

// Local stand-in for a serde-style deserializer entry point (borrows its input).
fn from_str(_input: &str) -> u32 {
    0
}

// FAILING: value cloned only to be deserialized.
fn wasteful(text: String) {
    let _ = from_str(&text.clone()); // should warn: value_cloned_to_deserialize
}

// PASSING: pass a reference to the original.
fn frugal(text: String) {
    let _ = from_str(&text);
}

fn main() {}
