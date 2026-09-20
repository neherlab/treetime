// The test dylint.toml bans `Vec::push` and `Option::expect` from `two_bans`, and `Vec::push`
// alone from `one_ban_twice`. Together these pin the counting rule: a finding is one
// (root, banned definition), not one root and not one call site.

// Two findings, one per banned definition. Before the walk ran to exhaustion this reported
// only the first ban it reached, leaving the other unwitnessed until that one was fixed.
fn two_bans(n: u32) -> u32 {
    let mut v = Vec::new();
    v.push(n);
    v[0] + Some(n).expect("n")
}

// One finding: the same definition reached twice is one ban broken, and the count of call
// sites rides on the finding rather than doubling it.
fn one_ban_twice(n: u32) -> u32 {
    let mut v = Vec::new();
    v.push(n);
    v.push(n + 1);
    v[0]
}

fn main() {
    let _ = two_bans(1);
    let _ = one_ban_twice(2);
}
