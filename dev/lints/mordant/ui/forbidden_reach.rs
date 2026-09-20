// The test dylint.toml bans Vec::push from anything matching `hot_path`.

// Flagged, with the witness path hot_path -> helper -> Vec::push.
fn hot_path(n: u32) -> u32 {
    helper(n)
}

fn helper(n: u32) -> u32 {
    let mut v = Vec::new();
    v.push(n);
    v[0]
}

// Fine: not a declared root, pushes all it likes.
fn cold_path(n: u32) -> u32 {
    let mut v = Vec::new();
    v.push(n);
    v[0]
}

// Fine: a root that never reaches the banned call.
fn hot_path_clean(n: u32) -> u32 {
    n + 1
}

// Flagged: `arr[i]` on a slice is a compiler builtin with no HIR call to see
// -- its bounds check exists only as a MIR `Assert` -- so the ban is broken
// here without ever passing through `callee_of`'s `Call`/`MethodCall` arms.
fn index_root(arr: &[u32], i: usize) -> u32 {
    arr[i]
}

fn main() {
    let _ = hot_path(1);
    let _ = cold_path(2);
    let _ = hot_path_clean(3);
    let _ = index_root(&[1, 2, 3], 0);
}
