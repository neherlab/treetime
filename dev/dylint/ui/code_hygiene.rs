#![allow(dead_code, unused_imports, unused_variables, clippy::all, reason = "ui fixture exercises patterns that are never executed")]

// FAILING: `use super::` import.
mod outer {
    pub struct Thing;
    pub mod inner {
        use super::Thing; // should warn: super_import
    }
}

// PASSING: absolute path from `crate::`.
mod outer_ok {
    pub struct Thing;
    pub mod inner {
        use crate::outer_ok::Thing; // no warning
    }
}

// FAILING: `use` inside a function; nested `fn`; both fire.
fn with_local_items() {
    use std::collections::BTreeMap; // should warn: local_use
    fn helper() {} // should warn: nested_function
    let _ = BTreeMap::<u8, u8>::new();
    helper();
}

// PASSING: closure capturing local state is allowed, not a nested fn.
fn with_closure(threshold: i32) -> i32 {
    let above = |x: i32| x > threshold;
    if above(3) { 1 } else { 0 }
}

// FAILING: version- and age-suffixed names.
fn parse_v2() {} // should warn: versioned_name
struct ConfigOld; // should warn: versioned_name
fn compute_new() {} // should warn: versioned_name

// PASSING: descriptive names, and a shadowed binding that is not an item name.
fn parse_header() {
    let value = 1;
    let value = value + 1; // shadowing is fine; not an item name
    let _ = value;
}

fn main() {}
