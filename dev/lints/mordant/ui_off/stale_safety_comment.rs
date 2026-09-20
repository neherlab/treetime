// The same shape ui/stale_safety_comment.rs flags; without
// `stale-safety-comment-enabled` nothing is reported.

fn stale_name(p: *const u64) -> u64 {
    // SAFETY: `frames_lock` is held for the duration of this read.
    unsafe { *p }
}

fn main() {
    let x = 7u64;
    println!("{}", stale_name(&x));
}
