// A crate fn folds a `Result`'s typed error into a bare `false`/`None`, and a call drops even that.

struct SysError(i32);

fn sys_write(buf: &mut Vec<u8>, bytes: &[u8]) -> Result<usize, SysError> {
    if buf.len() + bytes.len() > 64 {
        return Err(SysError(28));
    }
    buf.extend_from_slice(bytes);
    Ok(bytes.len())
}

// Collapses in a tail `match`: `Err(_) => false`.
fn write_pidfile(buf: &mut Vec<u8>) -> bool {
    match sys_write(buf, b"1") {
        Ok(_) => true,
        Err(_) => false,
    }
}

// Flagged: statement position drops the `false`.
fn pidfile_dropped(buf: &mut Vec<u8>) {
    write_pidfile(buf);
}

// Flagged: `let _ =` drops it just the same.
fn pidfile_discarded(buf: &mut Vec<u8>) {
    let _ = write_pidfile(buf);
}

// Fine: the `false` is read.
fn pidfile_checked(buf: &mut Vec<u8>) -> u8 {
    if write_pidfile(buf) { 1 } else { 0 }
}

// Collapses as the returned value: `.is_ok()`.
fn set_mode(buf: &mut Vec<u8>) -> bool {
    sys_write(buf, b"m").is_ok()
}

// Flagged.
fn mode_discarded(buf: &mut Vec<u8>) {
    let _ = set_mode(buf);
}

// Collapses in an early return: `if r.is_err() { return false }`.
fn sync_dir(buf: &mut Vec<u8>) -> bool {
    let r = sys_write(buf, b"s");
    if r.is_err() {
        return false;
    }
    buf.push(0);
    true
}

// Flagged.
fn sync_dropped(buf: &mut Vec<u8>) {
    sync_dir(buf);
}

// Collapses in a `let .. else`.
fn reserve(buf: &mut Vec<u8>) -> bool {
    let Ok(n) = sys_write(buf, b"r") else {
        return false;
    };
    n > 0
}

// Flagged.
fn reserve_discarded(buf: &mut Vec<u8>) {
    let _ = reserve(buf);
}

// Collapses into `None` through `.ok()?`.
fn open_slot(buf: &mut Vec<u8>) -> Option<usize> {
    let n = sys_write(buf, b"o").ok()?;
    Some(n + 1)
}

// Flagged: the `None` is dropped.
fn slot_dropped(buf: &mut Vec<u8>) {
    open_slot(buf);
}

// Fine: the `None` is read.
fn slot_checked(buf: &mut Vec<u8>) -> usize {
    if let Some(n) = open_slot(buf) { n } else { 0 }
}

// Fine: `?` on the `Option` passes the `None` on; and `slot_chained` itself
// collapses no `Result`, so dropping ITS `None` is not this lint's.
fn slot_chained(buf: &mut Vec<u8>) -> Option<usize> {
    let n = open_slot(buf)?;
    Some(n * 2)
}

fn chained_dropped(buf: &mut Vec<u8>) {
    slot_chained(buf);
}

// Collapses into `None` in a tail `if let .. else`.
fn probe_slot(buf: &mut Vec<u8>) -> Option<usize> {
    if let Ok(n) = sys_write(buf, b"p") {
        Some(n)
    } else {
        None
    }
}

// Flagged.
fn probe_dropped(buf: &mut Vec<u8>) {
    probe_slot(buf);
}

// Fine: the `Err` arm looks at the error before answering `false`.
fn write_logged(buf: &mut Vec<u8>, log: &mut Vec<i32>) -> bool {
    match sys_write(buf, b"l") {
        Ok(_) => true,
        Err(e) => {
            log.push(e.0);
            false
        }
    }
}

fn logged_dropped(buf: &mut Vec<u8>, log: &mut Vec<i32>) {
    write_logged(buf, log);
}

// Fine: an `Err` arm that names a kind asked the error a question, and its
// sibling saw the rest.
fn write_unless_full(buf: &mut Vec<u8>, log: &mut Vec<i32>) -> bool {
    match sys_write(buf, b"k") {
        Ok(_) => true,
        Err(SysError(28)) => false,
        Err(e) => {
            log.push(e.0);
            false
        }
    }
}

fn unless_full_dropped(buf: &mut Vec<u8>, log: &mut Vec<i32>) {
    write_unless_full(buf, log);
}

// Fine: the same, spelt as two `if let`s.
fn write_unless_full_if(buf: &mut Vec<u8>, log: &mut Vec<i32>) -> bool {
    let r = sys_write(buf, b"k");
    if let Err(SysError(28)) = r {
        return false;
    }
    if let Err(e) = r {
        log.push(e.0);
        return false;
    }
    true
}

fn unless_full_if_dropped(buf: &mut Vec<u8>, log: &mut Vec<i32>) {
    write_unless_full_if(buf, log);
}

// Collapses in a `let .. else` directly inside a `loop` body.
fn fill(buf: &mut Vec<u8>) -> bool {
    loop {
        let Ok(n) = sys_write(buf, b"f") else {
            return false;
        };
        if n > 0 {
            break;
        }
    }
    true
}

// Flagged.
fn fill_dropped(buf: &mut Vec<u8>) {
    fill(buf);
}

// Collapses a `&str` error, which carried a message.
fn parse_digit(s: &str) -> Result<u8, &'static str> {
    match s.as_bytes() {
        [b @ b'0'..=b'9'] => Ok(b - b'0'),
        [] => Err("empty"),
        _ => Err("not a digit"),
    }
}

fn is_digit(s: &str) -> bool {
    parse_digit(s).is_ok()
}

// Flagged.
fn digit_discarded() {
    let _ = is_digit("x");
}

// Fine: every call reads the `bool`.
fn try_mark(buf: &mut Vec<u8>) -> bool {
    sys_write(buf, b"c").is_ok()
}

fn mark_checked(buf: &mut Vec<u8>) -> u8 {
    if try_mark(buf) { 1 } else { 0 }
}

// Fine: `binary_search`'s `Err(idx)` is an answer, not a failure.
fn is_sorted_in(haystack: &[u32], needle: u32) -> bool {
    haystack.binary_search(&needle).is_ok()
}

fn sorted_discarded(haystack: &[u32]) {
    let _ = is_sorted_in(haystack, 3);
}

// Fine: a `()` error carries no kind to lose.
fn ping(up: bool) -> Result<(), ()> {
    if up { Ok(()) } else { Err(()) }
}

fn pinged(up: bool) -> bool {
    ping(up).is_ok()
}

fn ping_discarded() {
    let _ = pinged(true);
}

// Fine: a zero-sized error has one value, so `false` already says all it did.
struct Full;

fn alloc_slot(n: usize) -> Result<usize, Full> {
    if n > 8 { Err(Full) } else { Ok(n) }
}

fn try_alloc(n: usize) -> bool {
    alloc_slot(n).is_ok()
}

fn alloc_discarded() {
    let _ = try_alloc(9);
}

// Fine: a trait method's signature is the trait's, not the impl's.
trait Sink {
    fn put(&mut self, byte: u8) -> bool;
}

impl Sink for Vec<u8> {
    fn put(&mut self, byte: u8) -> bool {
        sys_write(self, &[byte]).is_ok()
    }
}

fn put_dropped(buf: &mut Vec<u8>) {
    buf.put(1);
}

// Fine: an exported C signature cannot return a `Result`.
fn checked_len(n: usize) -> Result<usize, SysError> {
    if n > 64 { Err(SysError(22)) } else { Ok(n) }
}

extern "C" fn exported_len_ok(n: usize) -> bool {
    checked_len(n).is_ok()
}

fn exported_discarded() {
    let _ = exported_len_ok(3);
}

// Fine: `unwrap_or(false)` on a `Result<bool, _>` folds "could not tell"
// into "no"; that shape is `defaulted_failure`'s to judge.
fn probe(buf: &mut Vec<u8>) -> Result<bool, SysError> {
    Ok(sys_write(buf, b"q")? > 0)
}

fn folded(buf: &mut Vec<u8>) -> bool {
    probe(buf).unwrap_or(false)
}

fn folded_discarded(buf: &mut Vec<u8>) {
    let _ = folded(buf);
}

// Fine: the collapse happens inside a closure, whose signature is its own.
fn with_closure(buf: &mut Vec<u8>) -> bool {
    let attempt = |b: &mut Vec<u8>| sys_write(b, b"x").is_ok();
    attempt(buf)
}

fn closure_dropped(buf: &mut Vec<u8>) {
    with_closure(buf);
}

fn main() {
    let mut buf = Vec::new();
    let mut log = Vec::new();
    pidfile_dropped(&mut buf);
    pidfile_discarded(&mut buf);
    let _ = pidfile_checked(&mut buf);
    mode_discarded(&mut buf);
    sync_dropped(&mut buf);
    reserve_discarded(&mut buf);
    slot_dropped(&mut buf);
    let _ = slot_checked(&mut buf);
    chained_dropped(&mut buf);
    probe_dropped(&mut buf);
    logged_dropped(&mut buf, &mut log);
    unless_full_dropped(&mut buf, &mut log);
    unless_full_if_dropped(&mut buf, &mut log);
    fill_dropped(&mut buf);
    digit_discarded();
    let _ = mark_checked(&mut buf);
    sorted_discarded(&[1, 2, 3]);
    ping_discarded();
    alloc_discarded();
    put_dropped(&mut buf);
    exported_discarded();
    folded_discarded(&mut buf);
    closure_dropped(&mut buf);
}
