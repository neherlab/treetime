//! An integer the function received is bounded by a branch somewhere in the
//! body and handed to something that turns it into memory somewhere else,
//! and that site is not dominated by any branch that tests it.
//!
//! The analysis is over MIR. A *seed* is an integer-typed place rooted at a
//! parameter other than `self` (`len`, `hdr.len`, `frame.header.count`, a
//! closure's capture), identified by its field path with dereferences
//! elided, so `hdr.len`, `(*hdr).len`, `let n = hdr.len` and a `&hdr.len`
//! read through later are one seed. `self` is the type's own state, and a
//! value the body computes itself -- a call's result or payload, `buf.len()`,
//! `n.min(cap)`, a loop variable -- is not input; bun's first run showed a
//! comparison of either is nearly always about something other than
//! trusting it. A parameter place the body writes or hands out as `&mut` is
//! not a seed either, since the analysis is flow-insensitive over its value;
//! a write through a copy of the parameter's pointer (`let q = p;
//! (*q).len = ..`) is a write to the same place.
//!
//! Reads are followed through single-assignment copies and borrows exactly,
//! and otherwise as marks. A cast or an aggregate (a `Range`, a tuple, a
//! `Some`) keeps the seed's [`Mark::Value`]; a borrow of it is a
//! [`Mark::Pointer`], which a read through a dereference turns back into the
//! `Value` and which sizes nothing itself (`from_raw_parts(&hdr.len, 1)`
//! views the field, it does not trust it); arithmetic turns either into a
//! [`Mark::Derived`] (`off + len` is a different quantity), as does storing
//! the value into a place some other store gives something else (`let mut q
//! = p; q += 1`, `n = n.min(cap)` on one path, `n = buf.len()` on another:
//! flow-insensitively, what such a place holds at a use is unknown); an
//! ordering comparison of the value or a derived quantity turns it into a
//! [`Mark::Bounded`] bool and an equality test into a [`Mark::Compared`] one,
//! followed the same way through `!` and `&&`; a call's result carries
//! `Compared` for whatever the arguments carried, by value or by reference
//! (`buf.get(i)?` judges `i`, so do `fits(&i)` and `i.cmp(&cap)`). A *check*
//! is a `SwitchInt` on a marked bool, hoisted out of the constant branch a
//! `debug_assert!` wraps its test in (`hoist_out_of_assertions`) so it counts
//! from where it is written. A seed qualifies only if some check bounds it:
//! an equality test, a judging call or the range test of a `match` arm
//! (`1..=8 =>` dispatches on the value like `0 =>` does) alone says nothing
//! about a size, but once the seed is bounded somewhere, any of them
//! dominating a use does cover it. A *sink* is a call from [`SIZE_SINKS`] or
//! [`PTR_SINKS`] given the seed's `Value`; safe indexing and slicing were
//! sinks in the first bun run and were removed (TRIAGE.md), since a panic on
//! a value the caller vouched for is what most of those turned out to be. A
//! sink is reported when the seed qualifies and no check dominates the
//! sink's block (rustc's forward dominators).

use std::collections::{HashMap, HashSet};

use rustc_abi::FieldIdx;
use rustc_hir::intravisit::{FnKind, Visitor, walk_pat};
use rustc_hir::{Body as HirBody, FnDecl, Pat, PatKind};
use rustc_index::IndexVec;
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::mir::{
    AggregateKind, BasicBlock, BinOp, Body, BorrowKind, Local, Operand, Place, ProjectionElem,
    RawPtrKind, Rvalue, StatementKind, TerminatorKind, VarDebugInfoContents,
};
use rustc_middle::ty::{self, TyCtxt};
use rustc_span::def_id::LocalDefId;
use rustc_span::symbol::kw;
use rustc_span::{Span, Symbol};

use crate::baseline::emit_with_note;
use crate::mir_flow::{local_name, mir_for};

rustc_session::declare_lint! {
    /// Flags an integer that comes from the caller (a parameter other than
    /// `self`, or a field reached from one) passed to `split_at`,
    /// `get_unchecked`, `with_capacity`, `reserve`, `set_len`,
    /// `from_raw_parts`, `copy_nonoverlapping` or a raw-pointer
    /// `add`/`offset` with no check on that path, although the same
    /// function does bound it (`<`, `<=`, `>`, `>=`, also via arithmetic on
    /// it) on another path: the author knew it needed a bound and this path
    /// skips it. Reported at the use, and the note points at the bound.
    ///
    /// Silent when every such use is dominated by a branch that tests the
    /// value (a clamp, an early return, an equality test, a `debug_assert!`,
    /// a `?` on a call it was passed to, by value or by reference), when no
    /// branch in the body bounds it at all (a check in the caller is
    /// invisible here, and a `match` arm's range is a dispatch, not a bound),
    /// when what reaches the use is arithmetic on the value or a copy the
    /// body re-assigns on some path rather than the value itself, when the
    /// value is `self`'s or one the body computed (a call's result or
    /// payload, `buf.len()`, a loop variable), when the place is written or
    /// `&mut`-borrowed anywhere in the body, directly or through a copy of
    /// its pointer, and on safe indexing and slicing, which panic rather than
    /// corrupt.
    ///
    /// Opt-in: runs when `dylint.toml` sets `unchecked-input-len-enabled`.
    /// The residual noise on bun (TRIAGE.md) is a value the caller vouches
    /// for that the function also uses as the limit of something else, which
    /// nothing inside the function separates from a real miss.
    pub UNCHECKED_INPUT_LEN,
    Warn,
    "received length turned into memory on a path that never checked it"
}

rustc_session::declare_lint_pass!(UncheckedInputLen => [UNCHECKED_INPUT_LEN]);

/// Calls that turn an integer argument into memory or a bound on it.
const SIZE_SINKS: &[&str] = &[
    "split_at",
    "split_at_mut",
    "split_at_unchecked",
    "split_at_mut_unchecked",
    "get_unchecked",
    "get_unchecked_mut",
    "with_capacity",
    "reserve",
    "reserve_exact",
    "set_len",
    "from_raw_parts",
    "from_raw_parts_mut",
    "copy_nonoverlapping",
    "copy_from_nonoverlapping",
    "copy_to_nonoverlapping",
];

/// Offsets on a raw pointer receiver. `add` on anything else is arithmetic.
const PTR_SINKS: &[&str] = &[
    "add",
    "sub",
    "offset",
    "byte_add",
    "byte_sub",
    "byte_offset",
];

/// One step of a place's field path. Derefs are not steps.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Step {
    Field(u32),
    Variant(u32),
}

/// A place with its dereferences elided: what the lint calls one value.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Key {
    local: Local,
    path: Vec<Step>,
}

/// A local assigned exactly once from another place: reads of it are reads
/// of that place.
enum Alias<'tcx> {
    /// `_l = copy p` / `_l = move p` / `_l = deref_copy p`: `_l.x` is `p.x`.
    Value(Place<'tcx>),
    /// `_l = &p` / `&raw p`: `(*_l).x` is `p.x`; `_l` itself is the borrow.
    Address(Place<'tcx>),
}

/// What a value carries about a seed, by seed id.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Mark {
    /// The seed's value, possibly cast or wrapped: usable as a size.
    Value(usize),
    /// A borrow or raw pointer to the seed's value: reading through it yields
    /// the `Value`, a call given it has looked at the seed, and it sizes
    /// nothing itself.
    Pointer(usize),
    /// Arithmetic on the seed, or a place that holds the seed on one path and
    /// something else on another: not the seed's size, but a comparison of
    /// it is still a check on the seed.
    Derived(usize),
    /// A bool that depends on an ordering comparison of the seed: the branch
    /// on it bounds the seed.
    Bounded(usize),
    /// A bool that depends on an equality test of the seed: the branch on it
    /// is a check, but not evidence that the seed needed bounding.
    Compared(usize),
}

type Marks = HashSet<Mark>;

struct Analysis<'a, 'tcx> {
    tcx: TyCtxt<'tcx>,
    body: &'a Body<'tcx>,
    aliases: HashMap<Local, Alias<'tcx>>,
    /// Single-assignment bools holding a constant: a branch on one is not a
    /// decision (`if cfg!(debug_assertions)`), and a check inside it counts
    /// from the branch itself.
    static_bools: HashSet<Local>,
    /// Parameter places the body writes or lends mutably; every seed they
    /// cover is dead. `self` and its fields are here from the start: the
    /// receiver's state is the type's own, not this call's input.
    killed: HashMap<Local, Vec<Vec<Step>>>,
    /// Places that hold a seed's value after one store and something else
    /// after another (`q += 1`, `n = n.min(cap)`, `n = buf.len()`): a copy of
    /// a seed stored into one of these is a different quantity, since the
    /// analysis cannot tell which store a use sees. Found by one propagation
    /// and applied by the next, until no new ones appear.
    accumulators: HashSet<Key>,
    found_accumulators: HashSet<Key>,
    /// Spans of the body's range patterns: the comparisons a `match` lowers
    /// them to are dispatches on the value, not bounds.
    range_patterns: Vec<Span>,
    seeds: Vec<Key>,
    seed_ids: HashMap<Key, usize>,
    /// Places holding values derived from seeds, by base local.
    taint: HashMap<Local, Vec<(Vec<Step>, Marks)>>,
}

fn overlap(a: &[Step], b: &[Step]) -> bool {
    a.starts_with(b) || b.starts_with(a)
}

/// The field path of a place, and whether every projection was one the path
/// represents. An index or subslice ends the path inexactly: what was read is
/// some part of the prefix.
fn key_of(place: Place<'_>) -> (Key, bool) {
    let mut key = Key {
        local: place.local,
        path: Vec::new(),
    };
    for elem in place.projection.iter() {
        match elem {
            ProjectionElem::Deref => {}
            ProjectionElem::Field(f, _) => key.path.push(Step::Field(f.as_u32())),
            ProjectionElem::Downcast(_, v) => key.path.push(Step::Variant(v.as_u32())),
            _ => return (key, false),
        }
    }
    (key, true)
}

fn mutably_lent<'tcx>(rvalue: &Rvalue<'tcx>) -> Option<Place<'tcx>> {
    match rvalue {
        Rvalue::Ref(_, BorrowKind::Mut { .. }, p) => Some(*p),
        Rvalue::RawPtr(RawPtrKind::Mut, p) => Some(*p),
        _ => None,
    }
}

impl<'a, 'tcx> Analysis<'a, 'tcx> {
    fn new(tcx: TyCtxt<'tcx>, body: &'a Body<'tcx>, range_patterns: Vec<Span>) -> Self {
        // Every place the body writes or lends out mutably.
        let mut targets: Vec<Place<'tcx>> = Vec::new();
        for data in body.basic_blocks.iter() {
            for stmt in &data.statements {
                match &stmt.kind {
                    StatementKind::Assign(a) => {
                        let (dest, rvalue) = &**a;
                        targets.push(*dest);
                        targets.extend(mutably_lent(rvalue));
                    }
                    StatementKind::SetDiscriminant { place, .. } => targets.push(**place),
                    _ => {}
                }
            }
            match &data.terminator().kind {
                TerminatorKind::Call { destination, .. } => targets.push(*destination),
                TerminatorKind::Yield { resume_arg, .. } => targets.push(*resume_arg),
                _ => {}
            }
        }
        // Writes to a local's own storage; a write through it (`(*q).len =
        // ..`) is a write to what it points at, which `collect_kills` charges
        // to the resolved place. Arguments carry one implicit write.
        let mut writes: IndexVec<Local, u32> = IndexVec::from_elem_n(0, body.local_decls.len());
        for l in body.args_iter() {
            writes[l] += 1;
        }
        for p in &targets {
            if !p.is_indirect() {
                writes[p.local] += 1;
            }
        }

        let mut aliases = HashMap::new();
        let mut static_bools = HashSet::new();
        for data in body.basic_blocks.iter() {
            for stmt in &data.statements {
                let StatementKind::Assign(a) = &stmt.kind else {
                    continue;
                };
                let (dest, rvalue) = &**a;
                if !dest.projection.is_empty() || writes[dest.local] != 1 {
                    continue;
                }
                let alias = match rvalue {
                    Rvalue::Use(op, _) => {
                        if op.place().is_none() {
                            static_bools.insert(dest.local);
                        }
                        op.place().map(Alias::Value)
                    }
                    // A reborrow is a bitwise copy of a reference-like ADT.
                    Rvalue::CopyForDeref(p) | Rvalue::Reborrow(_, _, p) => Some(Alias::Value(*p)),
                    Rvalue::Ref(_, _, p) | Rvalue::RawPtr(_, p) => Some(Alias::Address(*p)),
                    _ => None,
                };
                if let Some(alias) = alias {
                    aliases.insert(dest.local, alias);
                }
            }
        }

        let mut this = Analysis {
            tcx,
            body,
            aliases,
            static_bools,
            killed: HashMap::new(),
            accumulators: HashSet::new(),
            found_accumulators: HashSet::new(),
            range_patterns,
            seeds: Vec::new(),
            seed_ids: HashMap::new(),
            taint: HashMap::new(),
        };
        this.collect_kills(targets);
        this
    }

    /// Writes and mutable borrows whose target, once aliases are resolved,
    /// lies under a parameter; plus the receiver itself.
    fn collect_kills(&mut self, targets: Vec<Place<'tcx>>) {
        let receiver = self
            .body
            .var_debug_info
            .iter()
            .filter_map(|info| match info.value {
                VarDebugInfoContents::Place(p) if info.name == kw::SelfLower => Some(p),
                _ => None,
            });
        for target in receiver.chain(targets) {
            // A bare write to an alias local is the alias's own definition,
            // not a write to what it reads.
            if target.projection.is_empty() && self.aliases.contains_key(&target.local) {
                continue;
            }
            let (key, _) = key_of(self.resolve(target));
            if self.is_arg(key.local) {
                self.killed.entry(key.local).or_default().push(key.path);
            }
        }
    }

    fn resolve(&self, mut place: Place<'tcx>) -> Place<'tcx> {
        // Alias chains are acyclic (each link is a single definition read
        // after it), so this bounds work, not correctness.
        for _ in 0..16 {
            match self.aliases.get(&place.local) {
                Some(Alias::Value(base)) => {
                    place = base.project_deeper(place.projection, self.tcx);
                }
                Some(Alias::Address(base)) => match place.projection.split_first() {
                    Some((ProjectionElem::Deref, rest)) => {
                        place = base.project_deeper(rest, self.tcx);
                    }
                    _ => return place,
                },
                None => return place,
            }
        }
        place
    }

    fn is_arg(&self, local: Local) -> bool {
        (1..=self.body.arg_count).contains(&local.as_usize())
    }

    fn killed(&self, key: &Key) -> bool {
        self.killed
            .get(&key.local)
            .is_some_and(|paths| paths.iter().any(|p| overlap(p, &key.path)))
    }

    /// The seeds a read of `place` yields: the place itself, when it is one,
    /// plus whatever derived values were stored into or around it. Reading a
    /// borrow itself (`&len` handed to a call or put in a tuple) yields
    /// pointers to what it borrows; a read that ends in a dereference yields
    /// what the pointers stored there point at.
    fn marks_of_place(&mut self, place: Place<'tcx>) -> Marks {
        let resolved = self.resolve(place);
        if resolved.projection.is_empty()
            && let Some(Alias::Address(borrowed)) = self.aliases.get(&resolved.local)
        {
            let borrowed = *borrowed;
            return point_at(self.marks_of_place(borrowed));
        }
        let mut out = Marks::new();
        let (key, exact) = key_of(resolved);
        if exact
            && self.is_arg(key.local)
            && !self.killed(&key)
            && place.ty(&self.body.local_decls, self.tcx).ty.is_integral()
        {
            out.insert(Mark::Value(self.intern(key.clone())));
        }
        if let Some(entries) = self.taint.get(&key.local) {
            for (path, seeds) in entries {
                if overlap(path, &key.path) {
                    out.extend(seeds);
                }
            }
        }
        if let Some(ProjectionElem::Deref) = resolved.projection.last() {
            out = out
                .into_iter()
                .map(|m| match m {
                    Mark::Pointer(id) => Mark::Value(id),
                    Mark::Value(_) | Mark::Derived(_) | Mark::Bounded(_) | Mark::Compared(_) => m,
                })
                .collect();
        }
        out
    }

    fn marks_of_operand(&mut self, op: &Operand<'tcx>) -> Marks {
        match op.place() {
            Some(p) => self.marks_of_place(p),
            None => Marks::new(),
        }
    }

    fn intern(&mut self, key: Key) -> usize {
        if let Some(&id) = self.seed_ids.get(&key) {
            return id;
        }
        let id = self.seeds.len();
        self.seeds.push(key.clone());
        self.seed_ids.insert(key, id);
        id
    }

    /// Records `marks` as stored at `dest` plus `extra` steps below it.
    /// Returns whether anything new was learned.
    fn store(&mut self, dest: Place<'tcx>, extra: &[Step], marks: Marks) -> bool {
        let (mut key, _) = key_of(self.resolve(dest));
        key.path.extend_from_slice(extra);
        let accumulates = self.accumulators.contains(&key);
        let held = self
            .taint
            .get_mut(&key.local)
            .and_then(|entries| entries.iter_mut().find(|(p, _)| *p == key.path));
        let Some((_, held)) = held else {
            if marks.is_empty() {
                return false;
            }
            let marks = if accumulates {
                derive(marks, Mark::Derived)
            } else {
                marks
            };
            self.taint
                .entry(key.local)
                .or_default()
                .push((key.path, marks));
            return true;
        };
        // Some earlier store gave this place something other than what this
        // one does (or nothing at all), and one of them was a seed's value:
        // which of the two a use sees is not known here.
        if !accumulates && *held != marks && held.iter().chain(&marks).any(|m| sizes(*m)) {
            self.found_accumulators.insert(key);
        }
        let before = held.len();
        if accumulates {
            held.extend(derive(marks, Mark::Derived));
        } else {
            held.extend(marks);
        }
        held.len() != before
    }

    /// Flow-insensitive propagation of seeds through the body's assignments
    /// and calls, to a fixpoint (a loop can define a value from one assigned
    /// later), repeated from scratch while it turns up new accumulators.
    fn propagate(&mut self) {
        loop {
            self.taint.clear();
            self.found_accumulators.clear();
            self.propagate_once();
            if self.found_accumulators.is_empty() {
                return;
            }
            let found = std::mem::take(&mut self.found_accumulators);
            self.accumulators.extend(found);
        }
    }

    fn propagate_once(&mut self) {
        loop {
            let mut changed = false;
            for data in self.body.basic_blocks.iter() {
                for stmt in &data.statements {
                    let StatementKind::Assign(a) = &stmt.kind else {
                        continue;
                    };
                    let (dest, rvalue) = &**a;
                    // An alias's own definition: reads of it resolve to its
                    // source directly. A store through it is a store there.
                    if dest.projection.is_empty() && self.aliases.contains_key(&dest.local) {
                        continue;
                    }
                    changed |= self.propagate_assign(*dest, rvalue);
                }
                // What a call makes of its arguments is not their size, but a
                // branch on it (`buf.get(i)?`, `n.checked_add(4)?`) has judged
                // them, so it covers later uses without qualifying them.
                if let TerminatorKind::Call {
                    args, destination, ..
                } = &data.terminator().kind
                {
                    let mut marks = Marks::new();
                    for arg in args.iter() {
                        marks.extend(self.marks_of_operand(&arg.node));
                    }
                    changed |= self.store(*destination, &[], judged(marks));
                }
            }
            if !changed {
                return;
            }
        }
    }

    fn propagate_assign(&mut self, dest: Place<'tcx>, rvalue: &Rvalue<'tcx>) -> bool {
        match rvalue {
            Rvalue::Use(op, _)
            | Rvalue::Repeat(op, _)
            | Rvalue::Cast(_, op, _)
            | Rvalue::WrapUnsafeBinder(op, _) => {
                let marks = self.marks_of_operand(op);
                self.store(dest, &[], marks)
            }
            Rvalue::UnaryOp(_, op) => {
                let marks = derive(self.marks_of_operand(op), Mark::Derived);
                self.store(dest, &[], marks)
            }
            Rvalue::BinaryOp(op, ops) => {
                let mut marks = self.marks_of_operand(&ops.0);
                marks.extend(self.marks_of_operand(&ops.1));
                let marks = match op {
                    BinOp::Lt | BinOp::Le | BinOp::Gt | BinOp::Ge | BinOp::Cmp => {
                        derive(marks, Mark::Bounded)
                    }
                    BinOp::Eq | BinOp::Ne => derive(marks, Mark::Compared),
                    _ => derive(marks, Mark::Derived),
                };
                self.store(dest, &[], marks)
            }
            Rvalue::Discriminant(p) => {
                let marks = derive(self.marks_of_place(*p), Mark::Derived);
                self.store(dest, &[], marks)
            }
            Rvalue::Aggregate(kind, ops) => {
                let mut changed = false;
                for (i, op) in ops.iter_enumerated() {
                    let s = self.marks_of_operand(op);
                    let steps = self.aggregate_steps(kind, i);
                    changed |= self.store(dest, &steps, s);
                }
                changed
            }
            // A borrow that is not an alias: stored into a field, or into a
            // local some other path re-assigns.
            Rvalue::Ref(_, _, p) | Rvalue::RawPtr(_, p) => {
                let marks = point_at(self.marks_of_place(*p));
                self.store(dest, &[], marks)
            }
            Rvalue::Reborrow(..) | Rvalue::CopyForDeref(_) | Rvalue::ThreadLocalRef(_) => false,
        }
    }

    /// Where operand `i` of an aggregate lands, relative to the destination.
    fn aggregate_steps(&self, kind: &AggregateKind<'tcx>, i: FieldIdx) -> Vec<Step> {
        match kind {
            AggregateKind::Adt(did, vidx, _, _, active_field) => {
                let adt = self.tcx.adt_def(*did);
                let field = active_field.unwrap_or(i);
                if adt.is_enum() {
                    vec![Step::Variant(vidx.as_u32()), Step::Field(field.as_u32())]
                } else {
                    vec![Step::Field(field.as_u32())]
                }
            }
            AggregateKind::Tuple
            | AggregateKind::Closure(..)
            | AggregateKind::Coroutine(..)
            | AggregateKind::CoroutineClosure(..) => vec![Step::Field(i.as_u32())],
            // Elements are read through an index, which resolves to the
            // whole array; a raw pointer's parts are never read back.
            AggregateKind::Array(_) | AggregateKind::RawPtr(..) => Vec::new(),
        }
    }

    /// The branches that test each seed, keyed by seed, for the seeds at
    /// least one branch bounds. A seed only ever tested for equality is not
    /// in the map: nothing says it needed bounding. A branch on the integer
    /// itself (`match n`) is a dispatch, not a check, and is not in it
    /// either; the `<=` pair an arm's `1..=8` lowers to is the same dispatch,
    /// so it covers the arms like any test but bounds nothing.
    fn checks(&mut self) -> HashMap<usize, Checks> {
        let mut out: HashMap<usize, Checks> = HashMap::new();
        for (bb, data) in self.body.basic_blocks.iter_enumerated() {
            let term = data.terminator();
            let TerminatorKind::SwitchInt { discr, .. } = &term.kind else {
                continue;
            };
            let marks = self.marks_of_operand(discr);
            if marks.is_empty() {
                continue;
            }
            let span = term.source_info.span;
            let dispatches = self.range_patterns.iter().any(|p| p.contains(span));
            let from = self.hoist_out_of_assertions(bb, span);
            for mark in marks {
                let (id, bounds) = match mark {
                    Mark::Bounded(id) => (id, !dispatches),
                    Mark::Compared(id) => (id, false),
                    Mark::Value(_) | Mark::Pointer(_) | Mark::Derived(_) => continue,
                };
                out.entry(id).or_default().push((from, span, bounds));
            }
        }
        out.retain(|_, branches| branches.iter().any(|(_, _, bounds)| *bounds));
        out
    }

    /// `debug_assert!(n < cap)` is `if true { if !(n < cap) { panic } }`: the
    /// comparison's block is one arm of a branch on a constant and dominates
    /// nothing after the assertion, so the check counts from the outermost
    /// branch on a constant whose source encloses it. On the way up, a branch
    /// with a single non-diverging arm (the first half of `debug_assert!(a <
    /// n && b < n)`, an `assert!`, a `let .. else { panic }`) is stepped
    /// over: it turns nothing away, so it changes what the check dominates
    /// neither way, and counting the check from it would cover a use that
    /// sits between the two. Any other branch ends the climb.
    fn hoist_out_of_assertions(&self, check: BasicBlock, check_span: Span) -> BasicBlock {
        let dominators = self.body.basic_blocks.dominators();
        let mut from = check;
        let mut cur = check;
        while let Some(idom) = dominators.immediate_dominator(cur) {
            let term = self.body.basic_blocks[idom].terminator();
            if let TerminatorKind::SwitchInt { discr, targets } = &term.kind {
                let is_static = match discr.place() {
                    None => true,
                    Some(p) => p.projection.is_empty() && self.static_bools.contains(&p.local),
                };
                if is_static && term.source_info.span.source_callsite().contains(check_span) {
                    from = idom;
                } else {
                    // Exactly one way on: with `match` arms besides the
                    // panicking one, the check is only on its own arm's way.
                    let ways_on: HashSet<BasicBlock> = targets
                        .all_targets()
                        .iter()
                        .copied()
                        .filter(|t| !self.diverges(*t))
                        .collect();
                    if ways_on.len() != 1 {
                        break;
                    }
                }
            }
            cur = idom;
        }
        from
    }

    /// Control entering `bb` never comes back: a panic call, possibly behind
    /// the blocks that build its message.
    fn diverges(&self, mut bb: BasicBlock) -> bool {
        for _ in 0..8 {
            match &self.body.basic_blocks[bb].terminator().kind {
                TerminatorKind::Unreachable | TerminatorKind::UnwindTerminate(_) => return true,
                TerminatorKind::Call { target: None, .. } => return true,
                TerminatorKind::Call {
                    target: Some(next), ..
                }
                | TerminatorKind::Goto { target: next } => bb = *next,
                _ => return false,
            }
        }
        false
    }

    /// The seeds whose value is in `place`, paired with `place` when it is
    /// the seed itself.
    fn sized_by(&mut self, place: Place<'tcx>) -> Vec<(usize, Option<Place<'tcx>>)> {
        let marks = self.marks_of_place(place);
        let (key, _) = key_of(self.resolve(place));
        sized_by_values(marks)
            .into_iter()
            .map(|id| (id, (self.seeds[id] == key).then_some(place)))
            .collect()
    }

    fn sinks(&mut self) -> Vec<Sink<'tcx>> {
        let mut out = Vec::new();
        for (bb, data) in self.body.basic_blocks.iter_enumerated() {
            let TerminatorKind::Call {
                func,
                args,
                fn_span,
                ..
            } = &data.terminator().kind
            else {
                continue;
            };
            let Some((callee, _)) = func.const_fn_def() else {
                continue;
            };
            let Some(name) = self.tcx.opt_item_name(callee) else {
                continue;
            };
            let sized_by = if SIZE_SINKS.contains(&name.as_str()) {
                0..args.len()
            } else if PTR_SINKS.contains(&name.as_str())
                && args
                    .first()
                    .is_some_and(|a| a.node.ty(&self.body.local_decls, self.tcx).is_raw_ptr())
            {
                1..args.len()
            } else {
                continue;
            };
            let mut seeds: Vec<(usize, Option<Place<'tcx>>)> = Vec::new();
            for arg in &args[sized_by] {
                let Some(place) = arg.node.place() else {
                    continue;
                };
                for (id, via) in self.sized_by(place) {
                    match seeds.iter_mut().find(|(seen, _)| *seen == id) {
                        Some((_, known)) => {
                            if known.is_none() {
                                *known = via;
                            }
                        }
                        None => seeds.push((id, via)),
                    }
                }
            }
            seeds.sort_by_key(|(id, _)| *id);
            if !seeds.is_empty() {
                out.push(Sink {
                    block: bb,
                    span: *fn_span,
                    callee: name,
                    seeds,
                });
            }
        }
        out
    }

    /// The name the use spells the seed by: from the first user-named local
    /// on the copy chain the operand came through (`n` in `let n = hdr.len;
    /// ..(n)`, `q.len` in `let q = p; ..((*q).len)`), else from the parameter
    /// itself (`hdr.len`, a closure's captured `len`).
    fn name(&self, id: usize, via: Option<Place<'tcx>>) -> String {
        let key = &self.seeds[id];
        let mut cur = via;
        while let Some(place) = cur {
            if place.local == key.local {
                break;
            }
            if let Some(name) = local_name(self.body, place.local) {
                let (read, _) = key_of(place);
                return self.spell(name.to_string(), place.local, &read.path, 0);
            }
            cur = match self.aliases.get(&place.local) {
                Some(Alias::Value(base)) => Some(*base),
                Some(Alias::Address(_)) | None => None,
            };
        }
        // Debug entries rooted at the parameter's own local: the parameter,
        // a pattern binding in it, or a closure's capture. Deepest wins.
        let mut best: Option<(usize, Symbol)> = None;
        for info in &self.body.var_debug_info {
            let VarDebugInfoContents::Place(p) = info.value else {
                continue;
            };
            if info.composite.is_some() || p.local != key.local {
                continue;
            }
            let (k, exact) = key_of(p);
            if exact && key.path.starts_with(&k.path) && best.is_none_or(|(d, _)| k.path.len() > d)
            {
                best = Some((k.path.len(), info.name));
            }
        }
        let (named_depth, out) = match best {
            Some((depth, name)) => (depth, name.to_string()),
            None => (0, format!("_{}", key.local.as_usize())),
        };
        self.spell(out, key.local, &key.path, named_depth)
    }

    /// `out`, which names the first `named_depth` steps of `path` below
    /// `local`, extended with the field names of the remaining steps.
    fn spell(&self, mut out: String, local: Local, path: &[Step], named_depth: usize) -> String {
        let mut ty = self.body.local_decls[local].ty;
        let mut variant = None;
        for (depth, step) in path.iter().enumerate() {
            while let Some(inner) = ty.builtin_deref(true) {
                ty = inner;
            }
            match step {
                Step::Variant(v) => variant = Some(rustc_abi::VariantIdx::from_u32(*v)),
                Step::Field(f) => {
                    let idx = FieldIdx::from_u32(*f);
                    let (label, next) = match ty.kind() {
                        ty::Adt(adt, args) => {
                            let var = match (variant.take(), adt.is_enum()) {
                                (Some(v), true) => adt.variant(v),
                                (None, false) => adt.non_enum_variant(),
                                (Some(_), false) | (None, true) => return out,
                            };
                            let Some(field) = var.fields.get(idx) else {
                                return out;
                            };
                            (
                                field.name.to_string(),
                                field.ty(self.tcx, args).skip_normalization(),
                            )
                        }
                        ty::Tuple(tys) => match tys.get(idx.as_usize()) {
                            Some(t) => (f.to_string(), *t),
                            None => return out,
                        },
                        ty::Closure(_, args) => {
                            match args.as_closure().upvar_tys().get(idx.as_usize()) {
                                Some(t) => (f.to_string(), *t),
                                None => return out,
                            }
                        }
                        _ => return out,
                    };
                    if depth >= named_depth {
                        out.push('.');
                        out.push_str(&label);
                    }
                    ty = next;
                }
            }
        }
        out
    }
}

/// The branches testing one seed: block, span, and whether the test is an
/// ordering comparison rather than an equality.
type Checks = Vec<(BasicBlock, Span, bool)>;

struct Sink<'tcx> {
    block: BasicBlock,
    span: Span,
    callee: Symbol,
    /// The seeds whose value sizes it, each with the operand it arrived in
    /// when that operand is the seed itself (rather than a range built from
    /// it), for naming.
    seeds: Vec<(usize, Option<Place<'tcx>>)>,
}

/// What an operation on marked operands produces: a seed's value or a
/// quantity derived from it becomes `to` of that seed; arithmetic on or a
/// comparison of a pointer to it is neither its value nor a test of it; a
/// bool that already records a test passes through (`!(a < b)`, `x && y` as
/// `BitAnd`).
fn derive(marks: Marks, to: fn(usize) -> Mark) -> Marks {
    marks
        .into_iter()
        .map(|m| match m {
            Mark::Value(id) | Mark::Derived(id) => to(id),
            Mark::Pointer(id) => Mark::Derived(id),
            Mark::Bounded(_) | Mark::Compared(_) => m,
        })
        .collect()
}

/// What a call's result carries: it has looked at every seed it was given,
/// by value, derived, or by reference.
fn judged(marks: Marks) -> Marks {
    marks
        .into_iter()
        .map(|m| match m {
            Mark::Value(id) | Mark::Derived(id) | Mark::Pointer(id) => Mark::Compared(id),
            Mark::Bounded(_) | Mark::Compared(_) => m,
        })
        .collect()
}

/// What a borrow of a place carrying `marks` carries: a way to the seed's
/// value or a quantity derived from it. A borrow of a bool that recorded a
/// test, or of a call's result (`iter.next()` on a range built from the
/// seed), says nothing about the seed.
fn point_at(marks: Marks) -> Marks {
    marks
        .into_iter()
        .filter_map(|m| match m {
            Mark::Value(id) => Some(Mark::Pointer(id)),
            Mark::Pointer(_) | Mark::Derived(_) => Some(m),
            Mark::Bounded(_) | Mark::Compared(_) => None,
        })
        .collect()
}

/// A mark that can reach a sink as a seed's size, directly or read through.
fn sizes(mark: Mark) -> bool {
    matches!(mark, Mark::Value(_) | Mark::Pointer(_))
}

/// The seeds whose value itself is in `marks`, in a stable order.
fn sized_by_values(marks: Marks) -> Vec<usize> {
    let mut ids: Vec<usize> = marks
        .into_iter()
        .filter_map(|m| match m {
            Mark::Value(id) => Some(id),
            Mark::Pointer(_) | Mark::Derived(_) | Mark::Bounded(_) | Mark::Compared(_) => None,
        })
        .collect();
    ids.sort_unstable();
    ids
}

/// The spans of every range pattern in a body.
struct RangePatterns(Vec<Span>);

impl<'tcx> Visitor<'tcx> for RangePatterns {
    fn visit_pat(&mut self, pat: &'tcx Pat<'tcx>) {
        if let PatKind::Range(..) = pat.kind {
            self.0.push(pat.span);
        }
        walk_pat(self, pat);
    }
}

impl<'tcx> LateLintPass<'tcx> for UncheckedInputLen {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        _kind: FnKind<'tcx>,
        _decl: &'tcx FnDecl<'tcx>,
        hir_body: &'tcx HirBody<'tcx>,
        span: Span,
        def_id: LocalDefId,
    ) {
        if span.from_expansion() {
            return;
        }
        let Some(mir) = mir_for(cx.tcx, def_id) else {
            return;
        };
        let body: &Body<'tcx> = &mir;
        let mut range_patterns = RangePatterns(Vec::new());
        range_patterns.visit_body(hir_body);
        let mut an = Analysis::new(cx.tcx, body, range_patterns.0);
        an.propagate();
        let checks = an.checks();
        if checks.is_empty() {
            return;
        }
        let dominators = body.basic_blocks.dominators();
        let mut findings: Vec<(Span, String, Symbol, Span)> = Vec::new();
        for sink in an.sinks() {
            if !dominators.is_reachable(sink.block) {
                continue;
            }
            for (id, via) in sink.seeds {
                let Some(cs) = checks.get(&id) else {
                    continue;
                };
                // A sink is a call terminator and a check a switch terminator,
                // so the two are never one block and plain dominance is exact.
                if cs
                    .iter()
                    .any(|(c, ..)| dominators.dominates(*c, sink.block))
                {
                    continue;
                }
                let Some(bound) = cs
                    .iter()
                    .filter(|(_, _, bounds)| *bounds)
                    .map(|(_, s, _)| *s)
                    .min_by_key(|s| s.lo())
                else {
                    continue;
                };
                findings.push((sink.span, an.name(id, via), sink.callee, bound));
            }
        }
        findings.sort_by_key(|(span, name, ..)| (span.lo(), name.clone()));
        for (span, name, callee, bound) in findings {
            emit_with_note(
                cx,
                UNCHECKED_INPUT_LEN,
                span,
                format!(
                    "`{name}` comes from the caller and reaches `{callee}` here with no check on this path. The function does bound `{name}` on another path"
                ),
                bound,
                "the bound that does not cover this path",
                format!(
                    "check `{name}` on every path before this call, or once before the paths split"
                ),
            );
        }
    }
}
