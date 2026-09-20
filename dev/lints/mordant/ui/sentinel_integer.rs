// An integer field one function tests against a sentinel and another indexes
// with or offsets by unchecked: the sentinel reaches that use.

use std::collections::HashMap;

const INVALID_SLOT: u32 = u32::MAX;

struct Entry {
    slot: u32,
    parent: u32,
    depth: i32,
    width: u32,
    rank: u32,
}

struct Table {
    names: Vec<&'static str>,
    entries: Vec<Entry>,
    by_slot: HashMap<u32, usize>,
}

impl Table {
    // Fine: the comparison is the check.
    fn name(&self, e: &Entry) -> Option<&'static str> {
        if e.slot == INVALID_SLOT {
            None
        } else {
            Some(self.names[e.slot as usize])
        }
    }

    // Flagged: indexes with `slot` and nothing in this body tests it.
    fn rename(&mut self, e: &Entry, to: &'static str) {
        self.names[e.slot as usize] = to;
    }

    // Fine: bounds-tested through a local before the index.
    fn name_or_empty(&self, e: &Entry) -> &'static str {
        let i = e.slot as usize;
        if i >= self.names.len() {
            return "";
        }
        self.names[i]
    }

    // Fine: `assert_ne!` is the comparison, spelled as a macro.
    fn name_asserted(&self, e: &Entry) -> &'static str {
        assert_ne!(e.slot, INVALID_SLOT);
        self.names[e.slot as usize]
    }

    // Fine: `.get` followed by the early return is the bounds test.
    fn name_probed(&self, e: &Entry) -> &'static str {
        if self.names.get(e.slot as usize).is_none() {
            return "";
        }
        self.names[e.slot as usize]
    }

    // Fine: a keyed lookup answers for a key that is not there.
    fn forget(&mut self, e: &Entry) -> Option<usize> {
        let _ = self.by_slot[&e.slot];
        self.by_slot.remove(&e.slot)
    }

    // Flagged: the sum carries `slot`, whichever side it is on.
    fn name_after(&self, e: &Entry) -> &'static str {
        self.names[e.width as usize + e.slot as usize]
    }

    fn detach(&mut self, i: usize) {
        self.entries[i].parent = u32::MAX;
    }

    fn is_root(&self, e: &Entry) -> bool {
        e.parent == u32::MAX
    }

    // Flagged: `parent` is `u32::MAX` for a root, and here it offsets a
    // pointer with no test in sight.
    fn parent_ptr(&self, e: &Entry) -> *const Entry {
        unsafe { self.entries.as_ptr().add(e.parent as usize) }
    }

    // Flagged: wrapping arithmetic decides nothing about the sentinel.
    fn parent_ptr_wrapping(&self, e: &Entry) -> *const Entry {
        self.entries.as_ptr().wrapping_add(e.parent as usize)
    }

    // Fine: calls the helper that tests `parent`.
    fn parent_of(&self, e: &Entry) -> Option<&Entry> {
        if self.is_root(e) {
            return None;
        }
        Some(&self.entries[e.parent as usize])
    }

    // Fine: the arm that indexes is the one the sentinel cannot reach.
    fn parent_matched(&self, e: &Entry) -> Option<&Entry> {
        match e.parent {
            u32::MAX => None,
            _ => Some(&self.entries[e.parent as usize]),
        }
    }

    // Fine: only ever called from a function that tested `parent` first.
    fn parent_unchecked(&self, e: &Entry) -> &Entry {
        &self.entries[e.parent as usize]
    }

    fn grandparent(&self, e: &Entry) -> Option<&Entry> {
        if e.parent != u32::MAX {
            self.parent_of(self.parent_unchecked(e))
        } else {
            None
        }
    }

    fn forget_depth(&mut self, i: usize) {
        self.entries[i].depth = -1;
    }

    fn depth_known(&self, e: &Entry) -> bool {
        e.depth != -1
    }

    fn depth_or_zero(&self, e: &Entry) -> i32 {
        if self.depth_known(e) { e.depth } else { 0 }
    }

    // Flagged: `depth` is `-1` when unknown, and the sum indexes.
    fn below(&self, e: &Entry) -> &Entry {
        &self.entries[(e.depth + 1) as usize]
    }

    // Fine: `matches!` against the sentinel is the test.
    fn below_matched(&self, e: &Entry) -> Option<&Entry> {
        if matches!(e.depth, -1) {
            return None;
        }
        Some(&self.entries[(e.depth + 1) as usize])
    }

    // Fine: arithmetic alone meets no memory.
    fn deeper(&self, e: &Entry) -> i32 {
        e.depth + 1
    }

    // Fine: the clamp is the author deciding what an out-of-range value does.
    fn below_clamped(&self, e: &Entry) -> &Entry {
        let i = (e.depth + 1).max(0) as usize;
        &self.entries[i]
    }

    // Fine: `width` is set to `MAX` to mean unbounded and only ever ordered
    // against, never tested for equality: a bound, not a missing value.
    fn unbound(&mut self, i: usize) {
        self.entries[i].width = u32::MAX;
    }

    fn slack(&self, e: &Entry, used: u32) -> u32 {
        e.width - used
    }

    // Fine: `.get` answers for the sentinel itself.
    fn ranked(&self, e: &Entry) -> Option<&&'static str> {
        if e.rank == u32::MAX {
            return None;
        }
        self.names.get(e.rank as usize)
    }

    fn ranked_loose(&self, e: &Entry) -> Option<&&'static str> {
        self.names.get(e.rank as usize)
    }
}

fn main() {
    let mut t = Table {
        names: vec!["a", "b"],
        entries: vec![Entry {
            slot: 0,
            parent: INVALID_SLOT,
            depth: 0,
            width: 1,
            rank: 0,
        }],
        by_slot: HashMap::new(),
    };
    let e = Entry {
        slot: 1,
        parent: 0,
        depth: -1,
        width: u32::MAX,
        rank: u32::MAX,
    };
    let _ = t.name(&e);
    t.rename(&e, "c");
    let _ = t.name_or_empty(&e);
    let _ = t.name_asserted(&e);
    let _ = t.name_probed(&e);
    t.by_slot.insert(1, 0);
    let _ = t.forget(&e);
    let _ = t.name_after(&t.entries[0]);
    t.detach(0);
    let _ = t.parent_ptr(&e);
    let _ = t.parent_ptr_wrapping(&e);
    let _ = t.parent_of(&e).is_some();
    let _ = t.parent_matched(&e).is_some();
    let _ = t.grandparent(&e).is_some();
    t.forget_depth(0);
    let _ = t.depth_or_zero(&e);
    let _ = t.below(&t.entries[0]).slot + t.below_clamped(&e).slot;
    let _ = t.below_matched(&e).is_some();
    let _ = t.deeper(&e);
    t.unbound(0);
    let _ = t.slack(&e, 0);
    let _ = t.ranked(&e);
    let _ = t.ranked_loose(&e);
}
