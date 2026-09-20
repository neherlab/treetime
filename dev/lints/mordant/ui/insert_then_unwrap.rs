// A lookup-unwrap proven by an insert just above, with nothing in between
// that could touch the map or key, re-fetches a value the code already had.

use std::collections::HashMap;

fn helper() {}

fn adjacent_is_flagged() -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(1, 2);
    let v = m.get(&1).unwrap();
    *v
}

fn pure_statement_between_is_flagged() -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(7, 8);
    let offset = 1 + 2;
    let v = m.get(&7).unwrap();
    *v + offset
}

fn call_between_is_fine() -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(1, 2);
    helper();
    *m.get(&1).unwrap()
}

fn remove_between_is_fine() -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(1, 2);
    m.remove(&1);
    *m.get(&1).unwrap_or(&0)
}

fn different_key_is_fine() -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(1, 2);
    *m.get(&3).unwrap_or(&0)
}

fn rebound_key_is_fine(k: u32) -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(k + 1, 0);
    m.insert(k, 2);
    let k = k + 1;
    let v = m.get(&k).unwrap();
    *v
}

fn rebound_map_is_fine(k: u32, other: HashMap<u32, u32>) -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(k, 2);
    let m = other;
    let v = m.get(&k).unwrap();
    *v
}

fn rebound_in_tuple_is_fine(k: u32) -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(k + 1, 0);
    m.insert(k, 2);
    let (k, other) = (k + 1, 1);
    let v = m.get(&k).unwrap();
    *v + other
}

fn negated_key_is_fine(k: i32) -> i32 {
    let mut m: HashMap<i32, i32> = HashMap::new();
    m.insert(k, 0);
    m.insert(-k, 2);
    let v = m.get(&k).unwrap();
    *v
}

fn deref_key_is_flagged(k: &u32) -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(*k, 2);
    let v = m.get(k).unwrap();
    *v
}

fn lookup_in_the_rebinding_let_is_flagged(k: u32) -> u32 {
    let mut m: HashMap<u32, u32> = HashMap::new();
    m.insert(k, 2);
    let k = m.get(&k).unwrap();
    *k
}

fn main() {
    let _ = adjacent_is_flagged()
        + pure_statement_between_is_flagged()
        + call_between_is_fine()
        + remove_between_is_fine()
        + different_key_is_fine()
        + rebound_key_is_fine(1)
        + rebound_map_is_fine(1, HashMap::from([(1, 0)]))
        + rebound_in_tuple_is_fine(1)
        + deref_key_is_flagged(&1)
        + lookup_in_the_rebinding_let_is_flagged(1);
    let _ = negated_key_is_fine(1);
}
