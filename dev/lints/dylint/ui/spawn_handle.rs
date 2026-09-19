#![allow(dead_code, unused_variables, clippy::all, reason = "ui fixture exercises spawns that are never executed")]

// FAILING: the JoinHandle is dropped immediately.
fn drops_handle() {
    std::thread::spawn(|| {}); // should warn: spawn_handle_dropped
    let _ = std::thread::spawn(|| {}); // should warn: spawn_handle_dropped
}

// PASSING: the handle is bound and joined.
fn joins_handle() {
    let handle = std::thread::spawn(|| 1_u32);
    let _result = handle.join();
}

fn main() {}
