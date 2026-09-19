#![allow(dead_code, unused_variables, clippy::all, reason = "ui fixture exercises impls that are never used")]

use std::fmt;

struct Point {
    x: i32,
    y: i32,
}

// FAILING: hand-written Display and Debug impls.
impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y) // should warn: handwritten_fmt_impl (Display)
    }
}

impl fmt::Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Point") // should warn: handwritten_fmt_impl (Debug)
    }
}

// PASSING: derived Debug and a named rendering method.
#[derive(Debug)]
struct Vector {
    dx: i32,
    dy: i32,
}

impl Vector {
    fn render(&self) -> String {
        format!("({}, {})", self.dx, self.dy)
    }
}

fn main() {}
