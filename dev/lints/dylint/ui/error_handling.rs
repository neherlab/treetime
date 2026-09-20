#![allow(dead_code, unused_variables, clippy::all, reason = "ui fixture exercises error handling that is never executed")]

fn fallible() -> Result<u32, String> {
    Ok(1)
}

// FAILING: a Result bound to `_`, and a Result whose error is defaulted away.
fn discards() {
    let _ = fallible(); // should warn: discarded_result
    let _value = fallible().unwrap_or(0); // should warn: default_masks_error
    let _also = fallible().unwrap_or_default(); // should warn: default_masks_error
    let _lazy = fallible().unwrap_or_else(|_| 0); // should warn: default_masks_error
}

// PASSING: propagate, handle, or assert the error.
fn handles() -> Result<(), String> {
    let value = fallible()?; // propagated
    match fallible() {
        Ok(v) => println_stub(v),
        Err(e) => return Err(e),
    }
    // `unwrap_or_default` on an Option is fine: no error is discarded.
    let opt: Option<u32> = None;
    let _renamed = opt.unwrap_or_default(); // no warning: receiver is Option
    Ok(())
}

fn println_stub(_: u32) {}

fn main() {}
