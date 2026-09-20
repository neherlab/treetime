//! Finds which locals may hold a reference into another local, and whether such
//! a reference is used where the call to the separate function forbids it.

use rustc_data_structures::fx::FxHashSet;
use rustc_hir::LangItem;
use rustc_index::bit_set::DenseBitSet;
use rustc_middle::mir::visit::{MutatingUseContext, NonMutatingUseContext, PlaceContext, Visitor};
use rustc_middle::mir::{
    self, BasicBlock, Local, Location, Operand, Place, Rvalue, StatementKind, TerminatorKind,
};
use rustc_middle::ty::{self, GenericArgKind, Ty, TyCtxt};

use super::classify::add_locals_computed_from;
use super::region_io::{LocalFacts, block_live, pointer_copies_of};
use crate::mir_flow::reads_any;

/// Whether a value of this type can hold a reference or raw pointer to a local. True if it has a
/// lifetime, or a raw pointer in it or in a field (`Iov`, `NonNull`, `RawWaker`). Fields of a
/// foreign type with a `Drop` impl are skipped: it frees what they point to (`Vec<u8>`, `Rc`).
fn may_hold_address<'tcx>(
    tcx: TyCtxt<'tcx>,
    typing_env: ty::TypingEnv<'tcx>,
    ty: Ty<'tcx>,
) -> bool {
    if has_lifetime(ty) {
        return true;
    }
    let mut pending = vec![ty];
    let mut seen: FxHashSet<Ty<'tcx>> = FxHashSet::default();
    while let Some(ty) = pending.pop() {
        for arg in ty.walk() {
            let GenericArgKind::Type(inner) = arg.kind() else {
                continue;
            };
            if !seen.insert(inner) {
                continue;
            }
            if seen.len() > 4096 {
                return true;
            }
            match *inner.kind() {
                ty::RawPtr(..) => return true,
                ty::Adt(def, _) if def.has_dtor(tcx) && !def.did().is_local() => {}
                ty::Adt(def, args) => {
                    for field in def.all_fields() {
                        match tcx.try_normalize_erasing_regions(typing_env, field.ty(tcx, args)) {
                            Ok(ty) => pending.push(ty),
                            Err(_) => return true,
                        }
                    }
                }
                _ => {}
            }
        }
    }
    false
}

/// Whether a callee given a value of this type could store a reference
/// through it: `&mut Vec<&u8>`, `&Cell<Option<&T>>`, not `&[u8]` or `fmt::Arguments`.
fn may_take_address<'tcx>(
    tcx: TyCtxt<'tcx>,
    typing_env: ty::TypingEnv<'tcx>,
    ty: Ty<'tcx>,
) -> bool {
    // (type, whether data reached through it is writable: not after a `&`)
    let mut pending = vec![(ty, true)];
    let mut seen: FxHashSet<(Ty<'tcx>, bool)> = FxHashSet::default();
    let mut types_left = 4096u32;
    while let Some((ty, writable)) = pending.pop() {
        if !seen.insert((ty, writable)) {
            continue;
        }
        types_left -= 1;
        if types_left == 0 {
            return true;
        }
        match *ty.kind() {
            ty::RawPtr(pointee, _) => {
                if may_hold_address(tcx, typing_env, pointee) {
                    return true;
                }
            }
            ty::Ref(_, pointee, mutability) if writable && mutability.is_mut() => {
                if may_hold_address(tcx, typing_env, pointee) {
                    return true;
                }
            }
            ty::Ref(_, pointee, _) => {
                if may_hold_address(tcx, typing_env, pointee) {
                    if !pointee.is_freeze(tcx, typing_env) {
                        return true;
                    }
                    pending.push((pointee, false));
                }
            }
            ty::Adt(def, args) => {
                if let Some(boxed) = ty.boxed_ty() {
                    if !writable && !boxed.is_freeze(tcx, typing_env) {
                        return true;
                    }
                    pending.push((boxed, writable));
                } else if !tcx.is_lang_item(def.did(), LangItem::FormatArguments) {
                    for field in def.all_fields() {
                        match tcx.try_normalize_erasing_regions(typing_env, field.ty(tcx, args)) {
                            Ok(ty) => pending.push((ty, writable)),
                            Err(_) => return true,
                        }
                    }
                }
            }
            ty::Tuple(tys) => pending.extend(tys.iter().map(|ty| (ty, writable))),
            ty::Array(element, _) | ty::Slice(element) | ty::Pat(element, _) => {
                pending.push((element, writable));
            }
            ty::Closure(_, args) => {
                pending.extend(
                    args.as_closure()
                        .upvar_tys()
                        .iter()
                        .map(|ty| (ty, writable)),
                );
            }
            ty::Alias(..) => {
                match tcx.try_normalize_erasing_regions(typing_env, ty::Unnormalized::new_wip(ty)) {
                    Ok(normal) if normal != ty => pending.push((normal, writable)),
                    _ => return true,
                }
            }
            ty::Bool
            | ty::Char
            | ty::Int(_)
            | ty::Uint(_)
            | ty::Float(_)
            | ty::Str
            | ty::Never
            | ty::FnDef(..)
            | ty::FnPtr(..) => {}
            _ => return true,
        }
    }
    false
}

/// `Box` and `&mut`: moving one ends borrows through it. Not `&T` or `*mut T`.
pub(super) fn owns_pointee(ty: Ty<'_>) -> bool {
    !(ty.is_raw_ptr() || ty.ref_mutability() == Some(mir::Mutability::Not))
}

/// Whether the borrow checker records a borrow of `place` against
/// `place.local`: `&mut *m` yes, `&*r` on `r: &String` no.
fn borrow_of_local<'tcx>(tcx: TyCtxt<'tcx>, body: &mir::Body<'tcx>, place: Place<'tcx>) -> bool {
    place.iter_projections().all(|(base, elem)| {
        !matches!(elem, mir::ProjectionElem::Deref) || owns_pointee(base.ty(body, tcx).ty)
    })
}

/// Whether reading `place` copies out a pointer whose uses count as uses of
/// `place.local`. Only `Derefer` and `ElaborateBoxDerefs` write such reads.
pub(super) fn copies_borrow<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    place: Place<'tcx>,
) -> bool {
    if place.projection.is_empty() {
        return false;
    }
    if body.local_decls[place.local].ty.boxed_ty().is_some() {
        return true;
    }
    let read = place.ty(body, tcx).ty;
    place.is_indirect()
        && borrow_of_local(tcx, body, place)
        && (read.ref_mutability() == Some(mir::Mutability::Mut) || read.boxed_ty().is_some())
}

/// `may_hold_address`, or a `Box` (`Derefer` copies one from behind a `&`).
pub(super) fn can_hold_borrow<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    local: Local,
) -> bool {
    let ty = body.local_decls[local].ty;
    may_hold_address(tcx, body.typing_env(tcx), ty) || ty.boxed_ty().is_some()
}

pub(super) struct BorrowHolders {
    pub(super) any: DenseBitSet<Local>,
    pub(super) exclusive: DenseBitSet<Local>,
    /// Such a reference may be held somewhere that is not a local (`v.push(&x)`).
    pub(super) escaped: bool,
    pub(super) escaped_exclusive: bool,
}

impl BorrowHolders {
    fn escape(&mut self, exclusive: bool) {
        self.escaped = true;
        self.escaped_exclusive |= exclusive;
    }
}

/// Whether one operand holds a found reference (its `bool`) while another
/// could store a reference, so code that receives both can keep the first.
fn can_store_held_reference<T>(
    mut operands: impl Iterator<Item = (bool, T)> + Clone,
    storable: impl Fn(T) -> bool,
) -> bool {
    let holding = operands.clone().filter(|&(holds, _)| holds).count();
    operands.any(|(holds, operand)| holding > usize::from(holds) && storable(operand))
}

/// The strongest reference a piece of MIR holds into the locals of interest.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Held {
    Nothing,
    Shared,
    Exclusive,
}

impl Held {
    fn of(holds: bool, exclusive: bool) -> Self {
        match (holds, exclusive) {
            (_, true) => Held::Exclusive,
            (true, false) => Held::Shared,
            (false, false) => Held::Nothing,
        }
    }
    fn holds(self) -> bool {
        self != Held::Nothing
    }
    fn exclusive(self) -> bool {
        self == Held::Exclusive
    }
}

/// Answers whether one piece of MIR produces a reference into a local in
/// `of`, or reads a local in `found`, and whether exclusively.
struct ReadsBorrow<'a, 'mir, 'tcx> {
    tcx: TyCtxt<'tcx>,
    body: &'mir mir::Body<'tcx>,
    of: &'a DenseBitSet<Local>,
    found: &'a BorrowHolders,
    held: Held,
}

impl<'a, 'mir, 'tcx> ReadsBorrow<'a, 'mir, 'tcx> {
    fn of(
        tcx: TyCtxt<'tcx>,
        body: &'mir mir::Body<'tcx>,
        of: &'a DenseBitSet<Local>,
        found: &'a BorrowHolders,
        visit: impl FnOnce(&mut Self),
    ) -> Held {
        let mut visitor = ReadsBorrow {
            tcx,
            body,
            of,
            found,
            held: Held::Nothing,
        };
        visit(&mut visitor);
        visitor.held
    }
}

impl<'tcx> Visitor<'tcx> for ReadsBorrow<'_, '_, 'tcx> {
    fn visit_place(&mut self, place: &Place<'tcx>, context: PlaceContext, location: Location) {
        use MutatingUseContext as M;
        use NonMutatingUseContext as N;
        let exclusive_borrow =
            matches!(context, PlaceContext::MutatingUse(M::Borrow | M::RawBorrow));
        if self.of.contains(place.local) {
            let (holds, exclusive) = match context {
                _ if exclusive_borrow => {
                    let of_local = borrow_of_local(self.tcx, self.body, *place);
                    (of_local, of_local)
                }
                PlaceContext::NonMutatingUse(N::SharedBorrow | N::FakeBorrow | N::RawBorrow) => {
                    (borrow_of_local(self.tcx, self.body, *place), false)
                }
                PlaceContext::NonMutatingUse(N::Copy | N::Move | N::Inspect) => {
                    let of_local = copies_borrow(self.tcx, self.body, *place);
                    let unique = place.ty(self.body, self.tcx).ty.ref_mutability()
                        == Some(mir::Mutability::Mut);
                    (of_local, of_local && unique)
                }
                _ => (false, false),
            };
            self.held = self.held.max(Held::of(holds, exclusive));
        }
        if exclusive_borrow && place.is_indirect() && self.found.any.contains(place.local) {
            self.held = Held::Exclusive;
        }
        self.super_place(place, context, location);
    }

    fn visit_local(&mut self, local: Local, _: PlaceContext, _: Location) {
        let found = Held::of(
            self.found.any.contains(local),
            self.found.exclusive.contains(local),
        );
        self.held = self.held.max(found);
    }
}

/// The locals in `blocks` that may hold a reference into a local in `of`.
/// One too many loses a finding. One missed reports a change that fails.
pub(super) fn borrow_holders<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    blocks: &DenseBitSet<BasicBlock>,
    of: &DenseBitSet<Local>,
) -> BorrowHolders {
    let typing_env = body.typing_env(tcx);
    let can_hold =
        |place: &Place<'tcx>| !place.is_indirect() && can_hold_borrow(tcx, body, place.local);
    let storable = |operand: &Operand<'tcx>| {
        let moved_away = matches!(operand, Operand::Move(place)
            if place.projection.is_empty()
                && !may_hold_address(tcx, typing_env, body.local_decls[place.local].ty));
        !moved_away && may_take_address(tcx, typing_env, operand.ty(&body.local_decls, tcx))
    };
    let empty = DenseBitSet::new_empty(body.local_decls.len());
    let mut found = BorrowHolders {
        any: empty.clone(),
        exclusive: empty,
        escaped: false,
        escaped_exclusive: false,
    };
    loop {
        let mut added = false;
        for block in blocks.iter() {
            let data = &body.basic_blocks[block];
            for statement in &data.statements {
                match &statement.kind {
                    StatementKind::Assign(assign) => {
                        let (dest, rvalue) = &**assign;
                        let held = ReadsBorrow::of(tcx, body, of, &found, |p| {
                            p.visit_rvalue(rvalue, Location::START);
                        });
                        let (holds, exclusive) = (held.holds(), held.exclusive());
                        // `Both { r: &x, v: &mut v }` reaches a call as one argument.
                        let built = match rvalue {
                            Rvalue::Aggregate(_, operands) if holds && operands.len() > 1 => {
                                let by_operand: Vec<bool> = operands
                                    .iter()
                                    .map(|operand| {
                                        ReadsBorrow::of(tcx, body, of, &found, |p| {
                                            p.visit_operand(operand, Location::START);
                                        })
                                        .holds()
                                    })
                                    .collect();
                                can_store_held_reference(
                                    by_operand.iter().copied().zip(operands),
                                    storable,
                                )
                            }
                            _ => false,
                        };
                        let beside = !dest.projection.is_empty()
                            && !dest.is_indirect()
                            && (holds || found.any.contains(dest.local))
                            && may_take_address(tcx, typing_env, body.local_decls[dest.local].ty);
                        if built || beside {
                            found.escape(exclusive || found.exclusive.contains(dest.local));
                        }
                        if !holds {
                            continue;
                        }
                        if dest.is_indirect() {
                            found.escape(exclusive);
                        } else if can_hold(dest) {
                            added |= found.any.insert(dest.local);
                            if exclusive {
                                added |= found.exclusive.insert(dest.local);
                            }
                        }
                    }
                    StatementKind::Intrinsic(_) => {
                        let held = ReadsBorrow::of(tcx, body, of, &found, |p| {
                            p.visit_statement(statement, Location::START);
                        });
                        if held.holds() {
                            found.escape(held.exclusive());
                        }
                    }
                    _ => {}
                }
            }
            let Some(terminator) = &data.terminator else {
                continue;
            };
            match &terminator.kind {
                TerminatorKind::Call {
                    func,
                    args,
                    destination,
                    ..
                } => {
                    let reads = |operand: &Operand<'tcx>| {
                        ReadsBorrow::of(tcx, body, of, &found, |p| {
                            p.visit_operand(operand, Location::START);
                        })
                    };
                    let by_func = reads(func);
                    let by_arg: Vec<Held> = args.iter().map(|arg| reads(&arg.node)).collect();
                    let holds = by_arg.iter().any(|held| held.holds());
                    let exclusive = by_arg.iter().any(|held| held.exclusive());
                    if by_func.holds() {
                        found.escape(by_func.exclusive());
                    }
                    if !holds {
                        continue;
                    }
                    let arg_holds = by_arg.iter().map(|held| held.holds());
                    let can_store =
                        can_store_held_reference(arg_holds.zip(args), |arg| storable(&arg.node));
                    if can_store || destination.is_indirect() {
                        found.escape(exclusive);
                    }
                    if can_hold(destination) {
                        added |= found.any.insert(destination.local);
                        if exclusive {
                            added |= found.exclusive.insert(destination.local);
                        }
                    }
                }
                TerminatorKind::TailCall { .. }
                | TerminatorKind::InlineAsm { .. }
                | TerminatorKind::Yield { .. } => {
                    let held = ReadsBorrow::of(tcx, body, of, &found, |p| {
                        p.visit_terminator(terminator, Location::START);
                    });
                    if held.holds() {
                        found.escape(held.exclusive());
                    }
                }
                _ => {}
            }
        }
        if !added {
            break;
        }
    }
    found
}

fn has_lifetime(ty: Ty<'_>) -> bool {
    ty.walk()
        .any(|arg| matches!(arg.kind(), GenericArgKind::Lifetime(_)))
}

/// A `&T` where `T` holds no address: it borrows a parameter as shared only.
fn shared_ref_only<'tcx>(tcx: TyCtxt<'tcx>, typing_env: ty::TypingEnv<'tcx>, ty: Ty<'tcx>) -> bool {
    matches!(*ty.kind(), ty::Ref(_, pointee, mir::Mutability::Not)
        if !may_hold_address(tcx, typing_env, pointee))
}

pub(super) struct AfterExit<'a, 'tcx> {
    pub(super) tcx: TyCtxt<'tcx>,
    pub(super) body: &'a mir::Body<'tcx>,
    pub(super) locals: &'a LocalFacts,
    pub(super) part: &'a DenseBitSet<BasicBlock>,
    pub(super) entry: BasicBlock,
    pub(super) exit: BasicBlock,
}

/// Finds whether an item uses `local`, or a pointer copied out of it, in a
/// way a live borrow forbids: any use if exclusive, else a write, move or `&mut`.
struct ConflictingUse<'a> {
    local: Local,
    through: &'a DenseBitSet<Local>,
    exclusive: bool,
    found: bool,
}

impl<'tcx> Visitor<'tcx> for ConflictingUse<'_> {
    fn visit_place(&mut self, place: &Place<'tcx>, context: PlaceContext, location: Location) {
        if place.local == self.local || place.is_indirect() && self.through.contains(place.local) {
            self.found |= match context {
                PlaceContext::NonUse(_) => false,
                PlaceContext::MutatingUse(_)
                | PlaceContext::NonMutatingUse(NonMutatingUseContext::Move) => true,
                PlaceContext::NonMutatingUse(_) => self.exclusive,
            };
        }
        self.super_place(place, context, location);
    }

    fn visit_local(&mut self, local: Local, context: PlaceContext, _: Location) {
        if local == self.local && self.exclusive && context.is_use() {
            self.found = true;
        }
    }
}

impl AfterExit<'_, '_> {
    /// Whether the body after `exit` uses a parameter in a forbidden way
    /// while a result that may borrow it is live.
    pub(super) fn result_borrows_param(
        &self,
        exclusive_locals: &DenseBitSet<Local>,
        params: &DenseBitSet<Local>,
        returns: &DenseBitSet<Local>,
    ) -> bool {
        let body = self.body;
        let typing_env = body.typing_env(self.tcx);
        let borrowing_results: Vec<Local> = returns
            .iter()
            .filter(|&local| has_lifetime(body.local_decls[local].ty))
            .collect();
        if borrowing_results.is_empty() {
            return false;
        }
        let mut after: Option<(DenseBitSet<BasicBlock>, bool)> = None;
        let mut one_param = DenseBitSet::new_empty(body.local_decls.len());
        for param in params.iter() {
            one_param.clear();
            one_param.insert(param);
            let holders = borrow_holders(self.tcx, body, self.part, &one_param);
            for &result in &borrowing_results {
                if result == param || !holders.any.contains(result) {
                    continue;
                }
                let exclusive = exclusive_locals.contains(param)
                    || !shared_ref_only(self.tcx, typing_env, body.local_decls[result].ty);
                let (blocks, reenters) = after.get_or_insert_with(|| self.blocks_after_exit());
                if self.used_while_borrowed(blocks, *reenters, param, result, exclusive) {
                    return true;
                }
            }
        }
        false
    }

    /// `exit` and the blocks after it, and whether they enter the part again.
    fn blocks_after_exit(&self) -> (DenseBitSet<BasicBlock>, bool) {
        let body = self.body;
        let mut seen = DenseBitSet::new_empty(body.basic_blocks.len());
        let mut reenters = false;
        seen.insert(self.exit);
        let mut pending = vec![self.exit];
        while let Some(block) = pending.pop() {
            let Some(terminator) = &body.basic_blocks[block].terminator else {
                continue;
            };
            for next in terminator.successors() {
                if body.basic_blocks[next].is_cleanup {
                    continue;
                }
                if self.part.contains(next) {
                    reenters = true;
                } else if seen.insert(next) {
                    pending.push(next);
                }
            }
        }
        (seen, reenters)
    }

    /// Whether a block in `after` uses `param` in a forbidden way while
    /// `result`, or a local holding its borrow, is live or stored elsewhere.
    fn used_while_borrowed(
        &self,
        after: &DenseBitSet<BasicBlock>,
        reenters: bool,
        param: Local,
        result: Local,
        exclusive: bool,
    ) -> bool {
        let body = self.body;
        let holding = locals_holding(body, after, result);
        let (live_throughout, store_reenters) = self.blocks_after_untracked_store(after, &holding);
        let through = pointer_copies_of(self.tcx, body, after, param);
        let holding_live_at = |block: BasicBlock| {
            self.locals
                .live
                .get(&block)
                .is_some_and(|live| holding.iter().any(|local| live.contains(local)))
        };
        if exclusive && (store_reenters || reenters && holding_live_at(self.entry)) {
            return true;
        }
        for block in after.iter() {
            let live_throughout = live_throughout.contains(block);
            if !live_throughout && !holding_live_at(block) {
                continue;
            }
            let data = &body.basic_blocks[block];
            let live = (!live_throughout).then(|| live_before(body, self.locals, block, &holding));
            let holding_live_before = |index: usize| live.as_ref().is_none_or(|live| live[index]);
            let mut param_use = ConflictingUse {
                local: param,
                through: &through,
                exclusive,
                found: false,
            };
            for (index, statement) in data.statements.iter().enumerate() {
                if !holding_live_before(index) {
                    continue;
                }
                param_use.visit_statement(statement, Location::START);
                if param_use.found {
                    return true;
                }
            }
            if let Some(terminator) = &data.terminator
                && holding_live_before(data.statements.len())
            {
                param_use.visit_terminator(terminator, Location::START);
                if param_use.found {
                    return true;
                }
            }
        }
        false
    }

    /// The blocks of `after` reachable from a store that `escapes` finds, where
    /// the borrow counts as live throughout, and whether they enter the part.
    fn blocks_after_untracked_store(
        &self,
        after: &DenseBitSet<BasicBlock>,
        holding: &DenseBitSet<Local>,
    ) -> (DenseBitSet<BasicBlock>, bool) {
        let body = self.body;
        let mut live_throughout = DenseBitSet::new_empty(body.basic_blocks.len());
        let mut pending: Vec<BasicBlock> = after
            .iter()
            .filter(|&block| escapes(self.tcx, body, block, holding))
            .collect();
        for &block in &pending {
            live_throughout.insert(block);
        }
        let mut reenters = false;
        while let Some(block) = pending.pop() {
            let Some(terminator) = &body.basic_blocks[block].terminator else {
                continue;
            };
            for mut next in terminator.successors() {
                if self.part.contains(next) {
                    reenters = true;
                    next = self.exit;
                }
                if after.contains(next) && live_throughout.insert(next) {
                    pending.push(next);
                }
            }
        }
        (live_throughout, reenters)
    }
}

/// Whether `block` stores a borrow held by a local in `holding` somewhere
/// other than a local (`*slot = chunk`, `parts.push(chunk)`).
fn escapes<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    block: BasicBlock,
    holding: &DenseBitSet<Local>,
) -> bool {
    let data = &body.basic_blocks[block];
    let holds = |operand: &Operand<'tcx>| {
        has_lifetime(operand.ty(&body.local_decls, tcx))
            && reads_any(holding, |uses| uses.visit_operand(operand, Location::START))
    };
    for statement in &data.statements {
        let escaped = match &statement.kind {
            StatementKind::Assign(assign) => {
                let (dest, rvalue) = &**assign;
                dest.is_indirect()
                    && has_lifetime(dest.ty(body, tcx).ty)
                    && reads_any(holding, |uses| uses.visit_rvalue(rvalue, Location::START))
            }
            StatementKind::Intrinsic(_) => reads_any(holding, |uses| {
                uses.visit_statement(statement, Location::START)
            }),
            _ => false,
        };
        if escaped {
            return true;
        }
    }
    let Some(terminator) = &data.terminator else {
        return false;
    };
    match &terminator.kind {
        TerminatorKind::Call {
            func,
            args,
            destination,
            ..
        } => {
            if holds(func) {
                return true;
            }
            let by_arg: Vec<bool> = args.iter().map(|arg| holds(&arg.node)).collect();
            if !by_arg.contains(&true) {
                return false;
            }
            if destination.is_indirect() && has_lifetime(destination.ty(body, tcx).ty) {
                return true;
            }
            let typing_env = body.typing_env(tcx);
            can_store_held_reference(by_arg.iter().copied().zip(args), |arg| {
                may_take_address(tcx, typing_env, arg.node.ty(&body.local_decls, tcx))
            })
        }
        TerminatorKind::TailCall { .. }
        | TerminatorKind::InlineAsm { .. }
        | TerminatorKind::Yield { .. } => reads_any(holding, |uses| {
            uses.visit_terminator(terminator, Location::START)
        }),
        _ => false,
    }
}

/// `result` plus every local of a type with a lifetime assigned from one in the set.
fn locals_holding(
    body: &mir::Body<'_>,
    blocks: &DenseBitSet<BasicBlock>,
    result: Local,
) -> DenseBitSet<Local> {
    let mut holding = DenseBitSet::new_empty(body.local_decls.len());
    holding.insert(result);
    add_locals_computed_from(
        body,
        || blocks.iter(),
        &mut holding,
        |dest| !dest.is_indirect() && has_lifetime(body.local_decls[dest.local].ty),
    );
    holding
}

fn live_before(
    body: &mir::Body<'_>,
    locals: &LocalFacts,
    block: BasicBlock,
    holding: &DenseBitSet<Local>,
) -> Vec<bool> {
    let mut state = DenseBitSet::new_empty(body.local_decls.len());
    let mut live = vec![false; body.basic_blocks[block].statements.len() + 1];
    block_live(
        body,
        block,
        &locals.live,
        &locals.returned,
        &mut state,
        |index, state| live[index] = holding.iter().any(|local| state.contains(local)),
    );
    live
}
