//! Counts the statements in each piece of the body and marks those that use a
//! generic parameter. Also finds locals whose value is a constant in each copy.

use std::ops::ControlFlow;

use rustc_data_structures::fx::FxHashMap;
use rustc_hir::def_id::DefId;
use rustc_index::IndexVec;
use rustc_index::bit_set::DenseBitSet;
use rustc_middle::mir::visit::{NonMutatingUseContext, PlaceContext, Visitor};
use rustc_middle::mir::{
    self, BasicBlock, Local, Location, Operand, Place, Rvalue, Statement, StatementKind,
    Terminator, TerminatorKind,
};
use rustc_middle::ty::{
    self, Ty, TyCtxt, TypeSuperVisitable, TypeVisitable, TypeVisitableExt, TypeVisitor,
};
use rustc_span::{ExpnKind, Span, Spanned};

use crate::mir_flow::{FlowGraph, mir_for, reads_any};

/// For one basic block: items that do something at run time, those among
/// them that use a parameter ("dependent"), and those written by hand.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct BlockFacts {
    pub(super) counted: u32,
    pub(super) dependent: u32,
    pub(super) hand_written: u32,
}

/// Not produced by a macro. A desugaring of hand-written code counts.
pub(super) fn hand_written(span: Span) -> bool {
    let mut ctxt = span.ctxt();
    loop {
        if ctxt.is_root() {
            return true;
        }
        let expn = ctxt.outer_expn_data();
        match expn.kind {
            ExpnKind::Root => return true,
            ExpnKind::Macro(..) => return false,
            ExpnKind::Desugaring(_) | ExpnKind::AstPass(_) => ctxt = expn.call_site.ctxt(),
        }
    }
}

pub(super) fn counted_statement(statement: &Statement<'_>) -> bool {
    !matches!(
        statement.kind,
        StatementKind::StorageLive(_)
            | StatementKind::StorageDead(_)
            | StatementKind::Nop
            | StatementKind::FakeRead(_)
            | StatementKind::PlaceMention(_)
            | StatementKind::AscribeUserType(..)
            | StatementKind::Coverage(_)
            | StatementKind::ConstEvalCounter
            | StatementKind::BackwardIncompatibleDropHint { .. }
    )
}

pub(super) fn counted_terminator(terminator: &Terminator<'_>) -> bool {
    !matches!(
        terminator.kind,
        TerminatorKind::Goto { .. }
            | TerminatorKind::Return
            | TerminatorKind::Unreachable
            | TerminatorKind::UnwindResume
            | TerminatorKind::UnwindTerminate(_)
            | TerminatorKind::FalseEdge { .. }
            | TerminatorKind::FalseUnwind { .. }
    )
}

/// Tests whether a type, statement or terminator uses a type or const parameter. A closure's type
/// lists every parameter, but only its captures and signature (per use) and body (cached) count.
pub(super) struct ParamUses<'tcx> {
    tcx: TyCtxt<'tcx>,
    closure_bodies: FxHashMap<DefId, bool>,
}

impl<'tcx> TypeVisitor<TyCtxt<'tcx>> for ParamUses<'tcx> {
    type Result = ControlFlow<()>;

    fn visit_ty(&mut self, ty: Ty<'tcx>) -> ControlFlow<()> {
        if !ty.has_non_region_param() {
            return ControlFlow::Continue(());
        }
        match *ty.kind() {
            ty::Param(_) => ControlFlow::Break(()),
            ty::Closure(def, args) if !self.closure_uses_param(def, args) => {
                ControlFlow::Continue(())
            }
            _ => ty.super_visit_with(self),
        }
    }

    fn visit_const(&mut self, ct: ty::Const<'tcx>) -> ControlFlow<()> {
        if !ct.has_non_region_param() {
            return ControlFlow::Continue(());
        }
        match ct.kind() {
            ty::ConstKind::Param(_) => ControlFlow::Break(()),
            _ => ct.super_visit_with(self),
        }
    }
}

impl<'tcx> ParamUses<'tcx> {
    pub(super) fn new(tcx: TyCtxt<'tcx>) -> Self {
        Self {
            tcx,
            closure_bodies: FxHashMap::default(),
        }
    }

    pub(super) fn any(&mut self, value: &impl TypeVisitable<TyCtxt<'tcx>>) -> bool {
        value.has_non_region_param() && value.visit_with(self).is_break()
    }

    fn closure_uses_param(&mut self, def: DefId, args: ty::GenericArgsRef<'tcx>) -> bool {
        let parts = args.as_closure();
        self.any(&parts.tupled_upvars_ty())
            || self.any(&parts.sig_as_fn_ptr_ty())
            || self.closure_body_uses_param(def)
    }

    fn closure_body_uses_param(&mut self, def: DefId) -> bool {
        if let Some(&known) = self.closure_bodies.get(&def) {
            return known;
        }
        // The closure's own body names its type: assume "no" while reading it.
        self.closure_bodies.insert(def, false);
        let uses = match def.as_local().and_then(|def| mir_for(self.tcx, def)) {
            Some(body) => self.body_uses_param(&body),
            None => true,
        };
        self.closure_bodies.insert(def, uses);
        uses
    }

    fn body_uses_param(&mut self, body: &mir::Body<'tcx>) -> bool {
        body.local_decls.iter().any(|decl| self.any(&decl.ty))
            || body.basic_blocks.iter().any(|data| {
                data.statements
                    .iter()
                    .any(|s| self.statement_uses_param(body, s))
                    || data.terminator.as_ref().is_some_and(|t| self.any(t))
            })
    }

    fn statement_uses_param(
        &mut self,
        body: &mir::Body<'tcx>,
        statement: &Statement<'tcx>,
    ) -> bool {
        if let Some((place, Rvalue::Aggregate(kind, operands))) = statement.kind.as_assign()
            && let mir::AggregateKind::Closure(def, args) = **kind
        {
            return self.any(place)
                || self.any(&Ty::new_closure(self.tcx, def, args))
                || self.any(operands);
        }
        self.promoted_load_uses_param(body, statement)
            .unwrap_or_else(|| self.any(statement))
    }

    /// For `_n = const promoted[k]` of this body (which names every parameter): whether `k` does.
    fn promoted_load_uses_param(
        &mut self,
        body: &mir::Body<'tcx>,
        statement: &Statement<'tcx>,
    ) -> Option<bool> {
        let (_, Rvalue::Use(Operand::Constant(constant), _)) = statement.kind.as_assign()? else {
            return None;
        };
        let mir::Const::Unevaluated(unevaluated, ty) = constant.const_ else {
            return None;
        };
        let index = unevaluated.promoted?;
        let def = body.source.def_id();
        if unevaluated.def != def {
            return None;
        }
        // `promoted_mir` does not steal the body `mir_for` reads.
        let promoted = self.tcx.promoted_mir(def).get(index)?;
        Some(self.any(&ty) || self.body_uses_param(promoted))
    }
}

pub(super) fn locals_using_param<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
) -> DenseBitSet<Local> {
    let mut uses = ParamUses::new(tcx);
    let mut found = DenseBitSet::new_empty(body.local_decls.len());
    for (local, decl) in body.local_decls.iter_enumerated() {
        if uses.any(&decl.ty) {
            found.insert(local);
        }
    }
    found
}

/// Sets `found` if an item uses a place whose final type uses a parameter: not `(*_1).len: u32`.
pub(super) struct PlaceParams<'a, 'tcx> {
    body: &'a mir::Body<'tcx>,
    uses: ParamUses<'tcx>,
    found: bool,
}

impl<'tcx> Visitor<'tcx> for PlaceParams<'_, 'tcx> {
    fn visit_place(&mut self, place: &Place<'tcx>, _: PlaceContext, _: Location) {
        let ty = place.ty(&self.body.local_decls, self.uses.tcx).ty;
        if self.uses.any(&ty) {
            self.found = true;
        }
    }
}

impl<'a, 'tcx> PlaceParams<'a, 'tcx> {
    pub(super) fn new(tcx: TyCtxt<'tcx>, body: &'a mir::Body<'tcx>) -> Self {
        PlaceParams {
            body,
            uses: ParamUses::new(tcx),
            found: false,
        }
    }

    fn statement_params(&mut self, statement: &Statement<'tcx>) -> (bool, bool) {
        self.found = false;
        self.visit_statement(statement, Location::START);
        (
            self.uses.statement_uses_param(self.body, statement),
            self.found,
        )
    }

    fn terminator_params(&mut self, terminator: &Terminator<'tcx>) -> (bool, bool) {
        self.found = false;
        self.visit_terminator(terminator, Location::START);
        (self.uses.any(terminator), self.found)
    }

    pub(super) fn dependent_statement(&mut self, statement: &Statement<'tcx>) -> bool {
        let (names_param, typed_place) = self.statement_params(statement);
        names_param || typed_place
    }

    pub(super) fn dependent_terminator(&mut self, terminator: &Terminator<'tcx>) -> bool {
        let (names_param, typed_place) = self.terminator_params(terminator);
        names_param || typed_place
    }

    /// Names a parameter only in a constant, callee or cast: a compile-time constant in each copy.
    fn const_only(&mut self, rhs: AssignedValue<'_, 'tcx>) -> bool {
        let (names_param, typed_place) = match rhs {
            AssignedValue::Value(statement, _) => self.statement_params(statement),
            AssignedValue::Call(terminator, ..) => self.terminator_params(terminator),
        };
        names_param && !typed_place
    }
}

pub(super) fn classify<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
) -> (IndexVec<BasicBlock, BlockFacts>, usize) {
    let mut places = PlaceParams::new(tcx, body);
    let facts: IndexVec<BasicBlock, BlockFacts> = body
        .basic_blocks
        .iter()
        .map(|data| {
            if data.is_cleanup {
                return BlockFacts::default();
            }
            let mut facts = BlockFacts::default();
            let mut count = |dependent: bool, span: Span| {
                facts.counted += 1;
                if dependent {
                    facts.dependent += 1;
                }
                if hand_written(span) {
                    facts.hand_written += 1;
                }
            };
            for statement in &data.statements {
                if !counted_statement(statement) {
                    continue;
                }
                count(
                    places.dependent_statement(statement),
                    statement.source_info.span,
                );
            }
            if let Some(terminator) = &data.terminator
                && counted_terminator(terminator)
            {
                count(
                    places.dependent_terminator(terminator),
                    terminator.source_info.span,
                );
            }
            facts
        })
        .collect();
    let total = facts.iter().map(|block| block.counted as usize).sum();
    (facts, total)
}

/// Locals whose value is a compile-time constant in each copy
/// (`size_of::<T>()`, `width_of(TAG)`), and `blocks_under_const_branch`.
pub(super) fn per_copy_consts<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    flow: &FlowGraph,
) -> (DenseBitSet<Local>, DenseBitSet<BasicBlock>) {
    let mut places = PlaceParams::new(tcx, body);
    let mut consts = DenseBitSet::new_empty(body.local_decls.len());
    for (dest, rhs) in assignments(body, body.basic_blocks.indices()) {
        if let Some(local) = assigned_local(dest)
            && places.const_only(rhs)
        {
            consts.insert(local);
        }
    }
    loop {
        if !consts.is_empty() {
            add_locals_computed_from(
                body,
                || body.basic_blocks.indices(),
                &mut consts,
                |dest| assigned_local(dest).is_some(),
            );
        }
        let branch_blocks = blocks_under_const_branch(body, flow, &consts);
        if !consts_assigned_under_const_branch(body, &branch_blocks, &mut consts) {
            return (consts, branch_blocks);
        }
    }
}

/// Adds locals assigned a value that reads no local (`match TAG { 0 => 1, _ => 4 }`).
fn consts_assigned_under_const_branch(
    body: &mir::Body<'_>,
    blocks_under_const_branch: &DenseBitSet<BasicBlock>,
    per_copy_consts: &mut DenseBitSet<Local>,
) -> bool {
    let every = DenseBitSet::new_filled(body.local_decls.len());
    let mut added = false;
    for (dest, rhs) in assignments(body, blocks_under_const_branch.iter()) {
        if let Some(local) = assigned_local(dest)
            && !per_copy_consts.contains(local)
            && !matches!(
                rhs,
                AssignedValue::Value(_, Rvalue::Use(Operand::RuntimeChecks(_), _))
            )
            && !rhs.reads_any(&every)
        {
            added |= per_copy_consts.insert(local);
        }
    }
    added
}

/// Until nothing changes. Destination places are not searched.
pub(super) fn add_locals_computed_from<'tcx, I: Iterator<Item = BasicBlock>>(
    body: &mir::Body<'tcx>,
    blocks: impl Fn() -> I,
    set: &mut DenseBitSet<Local>,
    accept: impl Fn(Place<'tcx>) -> bool,
) {
    loop {
        let mut added = false;
        for (dest, rhs) in assignments(body, blocks()) {
            if !set.contains(dest.local) && accept(dest) && rhs.reads_any(set) {
                added |= set.insert(dest.local);
            }
        }
        if !added {
            return;
        }
    }
}

/// `None` for a store through a pointer or into an array element.
fn assigned_local(place: Place<'_>) -> Option<Local> {
    place
        .projection
        .iter()
        .all(|elem| {
            !matches!(
                elem,
                mir::ProjectionElem::Deref
                    | mir::ProjectionElem::Index(_)
                    | mir::ProjectionElem::ConstantIndex { .. }
                    | mir::ProjectionElem::Subslice { .. }
            )
        })
        .then_some(place.local)
}

#[derive(Clone, Copy)]
enum AssignedValue<'a, 'tcx> {
    Value(&'a Statement<'tcx>, &'a Rvalue<'tcx>),
    Call(
        &'a Terminator<'tcx>,
        &'a Operand<'tcx>,
        &'a [Spanned<Operand<'tcx>>],
    ),
}

impl<'tcx> AssignedValue<'_, 'tcx> {
    fn visit(self, visitor: &mut impl Visitor<'tcx>) {
        match self {
            AssignedValue::Value(_, rvalue) => visitor.visit_rvalue(rvalue, Location::START),
            AssignedValue::Call(_, func, args) => {
                visitor.visit_operand(func, Location::START);
                for arg in args {
                    visitor.visit_operand(&arg.node, Location::START);
                }
            }
        }
    }

    fn reads_any(self, of: &DenseBitSet<Local>) -> bool {
        reads_any(of, |uses| self.visit(uses))
    }
}

fn assignments<'a, 'tcx>(
    body: &'a mir::Body<'tcx>,
    blocks: impl Iterator<Item = BasicBlock>,
) -> impl Iterator<Item = (Place<'tcx>, AssignedValue<'a, 'tcx>)> {
    blocks
        .map(|block| &body.basic_blocks[block])
        .filter(|data| !data.is_cleanup)
        .flat_map(|data| {
            let assigns = data.statements.iter().filter_map(|statement| {
                let StatementKind::Assign(assign) = &statement.kind else {
                    return None;
                };
                Some((assign.0, AssignedValue::Value(statement, &assign.1)))
            });
            let call = data.terminator.as_ref().and_then(|terminator| {
                let TerminatorKind::Call {
                    func,
                    args,
                    destination,
                    ..
                } = &terminator.kind
                else {
                    return None;
                };
                Some((*destination, AssignedValue::Call(terminator, func, args)))
            });
            assigns.chain(call)
        })
}

/// The blocks that run only when a branch on a per-copy constant (`if
/// FLAG`) goes one way. Each copy keeps only its own side.
pub(super) fn blocks_under_const_branch(
    body: &mir::Body<'_>,
    flow: &FlowGraph,
    per_copy_consts: &DenseBitSet<Local>,
) -> DenseBitSet<BasicBlock> {
    let per_copy = |discr: &Operand<'_>| match discr {
        Operand::Constant(constant) => constant.const_.has_non_region_param(),
        Operand::Copy(place) | Operand::Move(place) => reads_any(per_copy_consts, |uses| {
            uses.visit_place(
                place,
                PlaceContext::NonMutatingUse(NonMutatingUseContext::Inspect),
                Location::START,
            );
        }),
        Operand::RuntimeChecks(_) => false,
    };
    let mut decided = DenseBitSet::new_empty(body.basic_blocks.len());
    let mut pending: Vec<BasicBlock> = body
        .basic_blocks
        .iter_enumerated()
        .filter(|&(block, data)| {
            flow.contains(block)
                && matches!(
                    &data.terminator,
                    Some(Terminator { kind: TerminatorKind::SwitchInt { discr, .. }, .. })
                        if per_copy(discr)
                )
        })
        .map(|(block, _)| block)
        .collect();
    let mut expanded = DenseBitSet::new_empty(body.basic_blocks.len());
    while let Some(branch) = pending.pop() {
        if !expanded.insert(branch) {
            continue;
        }
        for block in flow.decides(branch) {
            if decided.insert(block) && flow.succs(block).len() >= 2 {
                pending.push(block);
            }
        }
    }
    decided
}
