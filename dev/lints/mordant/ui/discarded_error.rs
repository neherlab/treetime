// Statement-position `.ok()` silently discards the error.

fn fallible() -> Result<u32, u32> {
    Err(1)
}

fn statement_ok_is_flagged() {
    fallible().ok();
}

fn let_underscore_is_fine() {
    let _ = fallible();
}

fn used_ok_is_fine() -> Option<u32> {
    fallible().ok()
}

fn main() {
    statement_ok_is_flagged();
    let_underscore_is_fine();
    let _ = used_ok_is_fine();
}
