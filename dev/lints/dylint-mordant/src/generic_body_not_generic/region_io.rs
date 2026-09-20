//! Finds what the chosen statements read and produce, and refuses them when the
//! call that would replace them could not compile or would run slower.

use std::collections::VecDeque;

use rustc_data_structures::fx::FxHashMap;
use rustc_data_structures::work_queue::WorkQueue;
use rustc_index::IndexVec;
use rustc_index::bit_set::DenseBitSet;
use rustc_middle::mir::visit::{MutatingUseContext, NonMutatingUseContext, PlaceContext, Visitor};
use rustc_middle::mir::{
    self, BasicBlock, Local, Location, Operand, Place, RETURN_PLACE, Rvalue, Statement,
    StatementKind, Terminator, TerminatorEdges, TerminatorKind,
};
use rustc_middle::ty::{self, Ty, TyCtxt};
use rustc_mir_dataflow::Analysis;
use rustc_mir_dataflow::impls::MaybeLiveLocals;

use super::borrows::{AfterExit, borrow_holders, can_hold_borrow, copies_borrow, owns_pointee};
use crate::mir_flow::{reaching, reads_any};

pub(super) struct LocalFacts {
    /// Liveness: the locals a later statement may still read, at the start
    /// of each non-cleanup block, ignoring unwind edges.
    pub(super) live: FxHashMap<BasicBlock, DenseBitSet<Local>>,
    pub(super) returned: DenseBitSet<Local>,
    storage: IndexVec<Local, Storage>,
    non_cleanup: DenseBitSet<BasicBlock>,
    per_copy_consts: DenseBitSet<Local>,
    /// Drop flags: hidden `bool`s the compiler adds to record whether a
    /// local still needs its drop. Found by shape, so a few others match.
    drop_flags: DenseBitSet<Local>,
    /// Each `switchInt` on a drop flag, with the local its `otherwise` edge drops.
    flag_guards: Vec<(Local, Option<Local>)>,
}

#[derive(Default)]
struct Storage {
    blocks: Vec<BasicBlock>,
    has_dead: bool,
}

impl LocalFacts {
    pub(super) fn new(body: &mir::Body<'_>, per_copy_consts: &DenseBitSet<Local>) -> Self {
        let mut storage = IndexVec::from_fn_n(|_| Storage::default(), body.local_decls.len());
        let mut non_cleanup = DenseBitSet::new_empty(body.basic_blocks.len());
        // `local_info` is gone at this phase: a user variable is one with debug info.
        let mut drop_flags = DenseBitSet::new_empty(body.local_decls.len());
        for (local, decl) in body.local_decls.iter_enumerated().skip(body.arg_count + 1) {
            if decl.ty.is_bool() {
                drop_flags.insert(local);
            }
        }
        for var in &body.var_debug_info {
            if let mir::VarDebugInfoContents::Place(place) = var.value {
                drop_flags.remove(place.local);
            }
        }
        for (block, data) in body.basic_blocks.iter_enumerated() {
            for statement in &data.statements {
                if let StatementKind::Assign(assign) = &statement.kind
                    && !assign.0.is_indirect()
                    && !matches!(assign.1, Rvalue::Use(Operand::Constant(_), _))
                {
                    drop_flags.remove(assign.0.local);
                }
            }
            if let Some(terminator) = &data.terminator
                && let TerminatorKind::Call { destination, .. } = &terminator.kind
                && !destination.is_indirect()
            {
                drop_flags.remove(destination.local);
            }
            if data.is_cleanup {
                continue;
            }
            non_cleanup.insert(block);
            for statement in &data.statements {
                match statement.kind {
                    StatementKind::StorageLive(local) => storage[local].blocks.push(block),
                    StatementKind::StorageDead(local) => {
                        storage[local].blocks.push(block);
                        storage[local].has_dead = true;
                    }
                    _ => {}
                }
            }
        }
        let mut flag_guards = Vec::new();
        for data in body.basic_blocks.iter() {
            if let Some(terminator) = &data.terminator
                && let TerminatorKind::SwitchInt { discr, targets } = &terminator.kind
                && let Some(place) = discr.place()
                && !place.is_indirect()
                && drop_flags.contains(place.local)
            {
                let guarded = match &body.basic_blocks[targets.otherwise()].terminator {
                    Some(Terminator {
                        kind: TerminatorKind::Drop { place, .. },
                        ..
                    }) if !place.is_indirect() => Some(place.local),
                    _ => None,
                };
                flag_guards.push((place.local, guarded));
            }
        }
        let mut returned = DenseBitSet::new_empty(body.local_decls.len());
        returned.insert(RETURN_PLACE);
        Self {
            live: live_in(body, &non_cleanup, &returned),
            returned,
            storage,
            non_cleanup,
            per_copy_consts: per_copy_consts.clone(),
            drop_flags,
            flag_guards,
        }
    }

    fn live_at(&self, block: BasicBlock) -> DenseBitSet<Local> {
        match self.live.get(&block) {
            Some(live) => live.clone(),
            None => DenseBitSet::new_empty(self.storage.len()),
        }
    }

    /// All the local's storage markers are in `blocks` and one is a
    /// `StorageDead`, so no borrow of it is usable outside `blocks`.
    fn storage_within(&self, local: Local, blocks: &DenseBitSet<BasicBlock>) -> bool {
        let storage = &self.storage[local];
        storage.has_dead && storage.blocks.iter().all(|&b| blocks.contains(b))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct InnerSignature {
    pub(super) params: Vec<Local>,
    pub(super) returns: Vec<Local>,
}

/// The signature of a separate fn holding `blocks`, or `None` when the call replacing them
/// would not compile or would run slower. `blocks` are left only toward `exit` (`None`: return).
pub(super) fn inner_signature<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    locals: &LocalFacts,
    blocks: &DenseBitSet<BasicBlock>,
    entry: BasicBlock,
    exit: Option<BasicBlock>,
) -> Option<InnerSignature> {
    debug_assert!(blocks.contains(entry) && exit.is_none_or(|x| !blocks.contains(x)));
    let uses = LocalUses::of(tcx, body, blocks);
    let live_out = match exit {
        Some(exit) => locals.live_at(exit),
        None => {
            let mut returned = DenseBitSet::new_empty(body.local_decls.len());
            returned.insert(RETURN_PLACE);
            returned
        }
    };
    // Moved out on some paths only: later code may only drop it, cheaply.
    if let Some(exit) = exit {
        let typing_env = body.typing_env(tcx);
        for local in uses.moved.iter() {
            if live_out.contains(local)
                && !uses.written.contains(local)
                && (body.local_decls[local]
                    .ty
                    .has_significant_drop(tcx, typing_env)
                    || read_after(body, exit, local))
            {
                return None;
            }
        }
    }
    let mut returns = live_out;
    returns.intersect(&uses.written);
    // A drop flag the part sets and later code tests may guard only such a local.
    for flag in returns.iter().filter(|&l| locals.drop_flags.contains(l)) {
        let mut tests = locals
            .flag_guards
            .iter()
            .filter(|&&(f, _)| f == flag)
            .peekable();
        if tests.peek().is_none()
            || tests.any(|&(_, guarded)| {
                guarded.is_none_or(|v| !uses.moved.contains(v) || uses.written.contains(v))
            })
        {
            return None;
        }
    }
    returns.subtract(&locals.drop_flags);
    // Live at the entry when liveness is computed over the part alone.
    let params = live_in(body, blocks, &returns)
        .remove(&entry)
        .unwrap_or_else(|| DenseBitSet::new_empty(body.local_decls.len()));
    if params
        .iter()
        .chain(returns.iter())
        .any(|local| !nameable(body.local_decls[local].ty))
    {
        return None;
    }
    // A drop flag live at the entry: the part could drop an uninitialized local.
    if params.iter().any(|local| locals.drop_flags.contains(local)) {
        return None;
    }
    if params
        .iter()
        .any(|local| locals.per_copy_consts.contains(local))
    {
        return None;
    }
    // The call moves or borrows each parameter: E0502 if a live local already does.
    let mut live = locals.live_at(entry);
    live.union(&params);
    if live.iter().any(|local| can_hold_borrow(tcx, body, local)) {
        let mut outside = reaching(body, entry);
        outside.intersect(&locals.non_cleanup);
        outside.subtract(blocks);
        let live_holds = |set: &DenseBitSet<Local>| live.iter().any(|local| set.contains(local));
        let mut exclusive_params = params.clone();
        exclusive_params.intersect(&uses.exclusive);
        if !exclusive_params.is_empty() {
            let holders = borrow_holders(tcx, body, &outside, &exclusive_params);
            if holders.escaped || live_holds(&holders.any) {
                return None;
            }
        }
        let mut shared = params.clone();
        shared.subtract(&uses.exclusive);
        if !shared.is_empty() {
            let holders = borrow_holders(tcx, body, &outside, &shared);
            if holders.escaped_exclusive || live_holds(&holders.exclusive) {
                return None;
            }
        }
    }
    // A reference into a local the separate fn owns must not outlive the call.
    let mut owned_addressed = uses.addressed.clone();
    for local in uses.addressed.iter() {
        let by_ref =
            params.contains(local) && !uses.moved.contains(local) && !returns.contains(local);
        if by_ref || locals.storage_within(local, blocks) {
            owned_addressed.remove(local);
        }
    }
    if !owned_addressed.is_empty() {
        let holders = borrow_holders(tcx, body, blocks, &owned_addressed);
        if holders.escaped || returns.iter().any(|local| holders.any.contains(local)) {
            return None;
        }
    }
    // A returned reference into a parameter borrows the whole parameter.
    if let Some(exit) = exit {
        let after = AfterExit {
            tcx,
            body,
            locals,
            part: blocks,
            entry,
            exit,
        };
        if after.result_borrows_param(&uses.exclusive, &params, &returns) {
            return None;
        }
    }
    Some(InnerSignature {
        params: params.iter().collect(),
        returns: returns.iter().collect(),
    })
}

/// What the blocks of a part do to each local they name.
struct LocalUses {
    /// Assigned, a call destination, mutably borrowed, or dropped. Not via deref.
    written: DenseBitSet<Local>,
    moved: DenseBitSet<Local>,
    /// Its own storage is borrowed (`&x`, `&mut x.f`), not its target (`&(*x).f`).
    addressed: DenseBitSet<Local>,
    /// Not passable as `&`: `written`, `moved`, and `through` that `owns_pointee`.
    exclusive: DenseBitSet<Local>,
    through: DenseBitSet<Local>,
    owning: DenseBitSet<Local>,
}

impl LocalUses {
    fn of<'tcx>(
        tcx: TyCtxt<'tcx>,
        body: &mir::Body<'tcx>,
        blocks: &DenseBitSet<BasicBlock>,
    ) -> Self {
        let empty = DenseBitSet::new_empty(body.local_decls.len());
        let mut owning = empty.clone();
        for (local, decl) in body.local_decls.iter_enumerated() {
            if owns_pointee(decl.ty) {
                owning.insert(local);
            }
        }
        let mut uses = LocalUses {
            written: empty.clone(),
            moved: empty.clone(),
            addressed: empty.clone(),
            exclusive: empty.clone(),
            through: empty,
            owning,
        };
        let mut copies: Vec<(Local, Local)> = Vec::new();
        for block in blocks.iter() {
            let data = &body.basic_blocks[block];
            for statement in &data.statements {
                uses.visit_statement(statement, Location::START);
                copies.extend(inserted_pointer_copy(tcx, body, statement));
            }
            if let Some(terminator) = &data.terminator
                && !matches!(terminator.kind, TerminatorKind::Return)
            {
                uses.visit_terminator(terminator, Location::START);
            }
        }
        loop {
            let mut added = false;
            for &(copy, of) in &copies {
                if uses.through.contains(copy) {
                    added |= uses.through.insert(of);
                }
            }
            if !added {
                break;
            }
        }
        for local in uses.through.iter() {
            if uses.owning.contains(local) {
                uses.exclusive.insert(local);
            }
        }
        uses.exclusive.union(&uses.written);
        uses.exclusive.union(&uses.moved);
        uses
    }
}

/// `(copy, original)` when the statement copies a pointer and the borrow
/// checker counts uses of the copy against the original (`copies_borrow`).
fn inserted_pointer_copy<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    statement: &Statement<'tcx>,
) -> Option<(Local, Local)> {
    let (dest, rvalue) = statement.kind.as_assign()?;
    if dest.is_indirect() {
        return None;
    }
    let (Rvalue::CopyForDeref(place)
    | Rvalue::Use(Operand::Copy(place), _)
    | Rvalue::Cast(_, Operand::Copy(place), _)) = rvalue
    else {
        return None;
    };
    copies_borrow(tcx, body, *place).then_some((dest.local, place.local))
}

pub(super) fn pointer_copies_of<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    blocks: &DenseBitSet<BasicBlock>,
    root: Local,
) -> DenseBitSet<Local> {
    let copies: Vec<(Local, Local)> = blocks
        .iter()
        .flat_map(|block| &body.basic_blocks[block].statements)
        .filter_map(|statement| inserted_pointer_copy(tcx, body, statement))
        .collect();
    let mut found = DenseBitSet::new_empty(body.local_decls.len());
    found.insert(root);
    loop {
        let mut added = false;
        for &(copy, of) in &copies {
            if found.contains(of) {
                added |= found.insert(copy);
            }
        }
        if !added {
            break;
        }
    }
    found
}

impl<'tcx> Visitor<'tcx> for LocalUses {
    fn visit_place(&mut self, place: &Place<'tcx>, context: PlaceContext, _: Location) {
        use MutatingUseContext as M;
        use NonMutatingUseContext as N;
        if !place.is_indirect() {
            let local = place.local;
            match context {
                PlaceContext::MutatingUse(
                    M::Store | M::SetDiscriminant | M::AsmOutput | M::Call | M::Yield | M::Retag,
                ) => {
                    self.written.insert(local);
                }
                PlaceContext::MutatingUse(M::Borrow | M::RawBorrow) => {
                    self.written.insert(local);
                    self.addressed.insert(local);
                }
                PlaceContext::MutatingUse(M::Drop) => {
                    self.written.insert(local);
                    self.moved.insert(local);
                }
                PlaceContext::NonMutatingUse(N::Move) => {
                    self.moved.insert(local);
                }
                PlaceContext::NonMutatingUse(N::SharedBorrow | N::FakeBorrow | N::RawBorrow) => {
                    self.addressed.insert(local);
                }
                PlaceContext::MutatingUse(M::Projection)
                | PlaceContext::NonMutatingUse(
                    N::Inspect | N::Copy | N::PlaceMention | N::Projection,
                )
                | PlaceContext::NonUse(_) => {}
            }
        } else if matches!(context, PlaceContext::MutatingUse(_)) {
            self.through.insert(place.local);
        }
    }
}

/// The locals live at the start of each block in `blocks`, along normal edges
/// inside `blocks` only. Where control leaves or returns, `boundary` is live.
fn live_in(
    body: &mir::Body<'_>,
    blocks: &DenseBitSet<BasicBlock>,
    boundary: &DenseBitSet<Local>,
) -> FxHashMap<BasicBlock, DenseBitSet<Local>> {
    let empty = DenseBitSet::new_empty(body.local_decls.len());
    let preds = body.basic_blocks.predecessors();
    let mut start: FxHashMap<BasicBlock, DenseBitSet<Local>> =
        blocks.iter().map(|block| (block, empty.clone())).collect();
    let mut queue: WorkQueue<BasicBlock> = WorkQueue::with_none(body.basic_blocks.len());
    let mut order: Vec<BasicBlock> = blocks.iter().collect();
    order.reverse();
    for block in order {
        queue.insert(block);
    }
    let mut state = empty;
    while let Some(block) = queue.pop() {
        block_live(body, block, &start, boundary, &mut state, |_, _| {});
        if let Some(known) = start.get_mut(&block)
            && known.union(&state)
        {
            for &pred in &preds[block] {
                if blocks.contains(pred) {
                    queue.insert(pred);
                }
            }
        }
    }
    start
}

/// Calls `each(i, state)` with the set of locals live just before statement `i`.
pub(super) fn block_live(
    body: &mir::Body<'_>,
    block: BasicBlock,
    live_at: &FxHashMap<BasicBlock, DenseBitSet<Local>>,
    boundary: &DenseBitSet<Local>,
    state: &mut DenseBitSet<Local>,
    mut each: impl FnMut(usize, &DenseBitSet<Local>),
) {
    live_after_statements(body, block, live_at, boundary, state);
    let data = &body.basic_blocks[block];
    each(data.statements.len(), state);
    for (index, statement) in data.statements.iter().enumerate().rev() {
        let at = Location {
            block,
            statement_index: index,
        };
        MaybeLiveLocals::transfer_function(state).visit_statement(statement, at);
        each(index, state);
    }
}

/// Sets `state` to the locals live after the last statement of `block`.
fn live_after_statements(
    body: &mir::Body<'_>,
    block: BasicBlock,
    live_at: &FxHashMap<BasicBlock, DenseBitSet<Local>>,
    boundary: &DenseBitSet<Local>,
    state: &mut DenseBitSet<Local>,
) {
    let data = &body.basic_blocks[block];
    state.clear();
    if let Some(terminator) = &data.terminator {
        let mut leaves = matches!(terminator.kind, TerminatorKind::Return);
        for next in terminator.successors() {
            let data = &body.basic_blocks[next];
            if data.is_cleanup || (data.terminator.is_some() && data.is_empty_unreachable()) {
                continue;
            }
            match live_at.get(&next) {
                Some(live) => {
                    state.union(live);
                }
                None => leaves = true,
            }
        }
        if leaves {
            state.union(boundary);
        }
        if !matches!(terminator.kind, TerminatorKind::Return) {
            // Only the edge where a call comes back is followed: write first.
            if let TerminatorEdges::AssignOnReturn { return_, place, .. } = terminator.edges()
                && !return_.is_empty()
            {
                MaybeLiveLocals.apply_call_return_effect(state, block, place);
            }
            let at = Location {
                block,
                statement_index: data.statements.len(),
            };
            MaybeLiveLocals::transfer_function(state).visit_terminator(terminator, at);
        }
    }
}

/// Whether `from` or a block after it reads `local`, other than its `Drop`
/// and the `Discriminant` read that starts an enum's drop.
fn read_after(body: &mir::Body<'_>, from: BasicBlock, local: Local) -> bool {
    let mut of = DenseBitSet::new_empty(body.local_decls.len());
    of.insert(local);
    let mut seen = DenseBitSet::new_empty(body.basic_blocks.len());
    let mut queue = VecDeque::from([from]);
    seen.insert(from);
    while let Some(block) = queue.pop_front() {
        let data = &body.basic_blocks[block];
        let found = reads_any(&of, |uses| {
            for statement in &data.statements {
                if let StatementKind::Assign(assign) = &statement.kind
                    && let Rvalue::Discriminant(place) = &assign.1
                    && !place.is_indirect()
                    && place.local == local
                {
                    continue;
                }
                uses.visit_statement(statement, Location::START);
            }
            if let Some(terminator) = &data.terminator
                && !matches!(
                    &terminator.kind,
                    TerminatorKind::Drop { place, .. } if !place.is_indirect() && place.local == local
                )
            {
                uses.visit_terminator(terminator, Location::START);
            }
        });
        if found {
            return true;
        }
        if let Some(terminator) = &data.terminator {
            for next in terminator.successors() {
                if !body.basic_blocks[next].is_cleanup && seen.insert(next) {
                    queue.push_back(next);
                }
            }
        }
    }
    false
}

fn nameable(ty: Ty<'_>) -> bool {
    !ty.walk().any(|arg| {
        matches!(
            arg.as_type().map(|ty| ty.kind()),
            Some(
                ty::Closure(..)
                    | ty::CoroutineClosure(..)
                    | ty::Coroutine(..)
                    | ty::CoroutineWitness(..)
                    | ty::FnDef(..)
                    | ty::Alias(ty::AliasTy {
                        kind: ty::Opaque { .. },
                        ..
                    })
            )
        )
    })
}
