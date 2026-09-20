//! Dataflow over a constructor's MIR: which fields of the returned `Self`
//! does a failure exit actually depend on?
//!
//! A receiver-less `fn(..) -> Result<Self, E>` / `Option<Self>` can fail for
//! three unrelated reasons: it inspected a value it is about to store and
//! rejected it (an invariant), the input it was reading was malformed (a
//! parser), or the environment refused (allocation, IO). Only the first makes
//! a field's later reassignment dangerous, and the signature cannot tell them
//! apart. The body can:
//!
//! * every point that returns `Err`/`None`, literal or via `?`, is a failure
//!   exit;
//! * the branch conditions an exit is control-dependent on (post-dominator
//!   analysis, so nesting and loops are exact) decide it, and everything those
//!   conditions read, followed back through copies, projections, arithmetic
//!   and predicate calls into their arguments, is the exit's *decision slice*;
//! * every field of the `Self` that reaches the return has a *storage slice*:
//!   the values it is built from, followed back through copies, casts,
//!   arithmetic and nested aggregates, but not through "the `Ok` payload of"
//!   or "the result of a call" -- past those the field's own type or the
//!   callee carries whatever guarantee exists, not this function's check;
//! * a field carries an invariant iff some exit's decision slice meets its
//!   storage slice.
//!
//! Slices are over *places* (`opts.footer`, not `opts`), so a check on one
//! field of an input struct is not a check on its siblings, and `args.port =
//! x; if args.port == 0 { Err }; Ok(args)` implicates `port` alone.
//!
//! `let v = input.next()?; Self { v }` fails on `input` and stores the payload
//! of the checked thing, so the slices never meet. `if n > MAX { Err }; Self {
//! n }` reads and stores `n`. Exits whose error type is a resource error
//! (`AllocError`, `io::Error`, ...) are dropped up front, since "the allocator
//! said no" and "the helper said no" have the same shape.
//!
//! The parts of this that know nothing about constructors -- which MIR to
//! read, places as atoms, the pruned CFG and control dependence -- live in
//! `mir_flow`; this module owns the def facts, the two slices, the resource
//! error list, helper summaries and the per-field verdict.

use std::collections::{HashMap, HashSet, VecDeque};

use rustc_abi::FieldIdx;
use rustc_hir::LangItem;
use rustc_hir::def_id::{DefId, LocalDefId};
use rustc_index::IndexVec;
use rustc_lint::LateContext;
use rustc_middle::mir::{
    AggregateKind, BasicBlock, Body, BorrowKind, Local, Mutability, Operand, Place, RETURN_PLACE,
    RawPtrKind, Rvalue, StatementKind, TerminatorKind,
};
use rustc_middle::ty::{self, Ty, TyCtxt};
use rustc_span::{Span, sym};

use crate::adt_facts::{matches_config_path, result_err_ty};
use crate::mir_flow::{
    ANY_ELEM, Atom, Exactness, FlowGraph, mir_for, place_info, switch_operand_atoms,
};

/// Error types that report the environment refusing, not the value being
/// wrong. An exit failing with one of these is never a validation.
const RESOURCE_ERRORS: &[&str] = &[
    "core::alloc::AllocError",
    "std::alloc::AllocError",
    "alloc::collections::TryReserveError",
    "std::collections::TryReserveError",
    "core::alloc::LayoutError",
    "std::alloc::LayoutError",
    "std::io::Error",
    "std::thread::AccessError",
];

pub(crate) struct FieldCheck {
    pub field: FieldIdx,
    /// The branch whose outcome the field's stored value decides.
    pub check: Span,
}

/// The fields of `self_did` that `ctor` checks before constructing one.
pub(crate) fn checked_fields(
    cx: &LateContext<'_>,
    ctor: LocalDefId,
    self_did: DefId,
    extra_resource_errors: &[String],
) -> Vec<FieldCheck> {
    let mut memo = HashMap::new();
    let mut out: HashMap<FieldIdx, Span> = HashMap::new();
    analyze_body(
        cx.tcx,
        ctor,
        self_did,
        extra_resource_errors,
        &mut memo,
        0,
        &mut out,
    );
    let mut v: Vec<FieldCheck> = out
        .into_iter()
        .map(|(field, check)| FieldCheck { field, check })
        .collect();
    v.sort_by_key(|f| f.field);
    v
}

/// The branch in `callee` whose outcome sends it to a non-resource
/// `Err`/`None` and whose condition reads back to one of its arguments other
/// than a `self` receiver: the check a caller discards when it replaces the
/// failure with a default. Only the branches an exit is directly
/// control-dependent on are read: an argument test that merely stands in
/// front of a receiver-decided exit (`if s.is_empty() { return Ok(0) }`
/// before `if self.closed { return Err(..) }`) accepts the argument, it does
/// not reject it. None when the body always succeeds, fails only on
/// resources, on its receiver's own state (an empty queue, a lexer at the
/// wrong token) or on state it was not handed at all, or builds its failure
/// with combinators instead of a branch of its own.
pub(crate) fn argument_decided_failure(
    tcx: TyCtxt<'_>,
    callee: LocalDefId,
    extra_resource_errors: &[String],
) -> Option<Span> {
    let first_argument = if tcx
        .opt_associated_item(callee.to_def_id())
        .is_some_and(|item| item.is_method())
    {
        2
    } else {
        1
    };
    let body = mir_for(tcx, callee)?;
    if body.tainted_by_errors.is_some() {
        return None;
    }
    if let Some(e) = result_err_ty(tcx, body.local_decls[RETURN_PLACE].ty)
        && is_resource_error(tcx, e, extra_resource_errors)
    {
        return None;
    }
    let facts = gather(tcx, &body, None, extra_resource_errors);
    if facts.failure_blocks.is_empty() {
        return None;
    }
    let flow = FlowGraph::new(&body);
    let mut visited: HashSet<BasicBlock> = HashSet::new();
    for &fb in &facts.failure_blocks {
        for branch in flow.direct_control_deps(fb) {
            if !visited.insert(branch)
                || decides_on_resource(tcx, &body, &facts, branch, extra_resource_errors)
            {
                continue;
            }
            let decision = decision_slice(&facts.defs, &switch_operand_atoms(&body, branch), &[]);
            if decision
                .iter()
                .any(|a| (first_argument..=body.arg_count).contains(&a.local.as_usize()))
            {
                return Some(body.basic_blocks[branch].terminator().source_info.span);
            }
        }
    }
    None
}

// ── def facts ────────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
enum Read {
    /// Value identity: a copy/move/borrow of the place. Paths compose.
    Same(Atom),
    /// Computed from the whole operand (cast, arithmetic, repeat).
    Derived(Atom),
    /// Operand `i` of an aggregate written to the dest. Paths compose past `i`.
    AggField(u32, Atom),
    /// A variant payload of this atom was read (`(x as Some).0`).
    Payload(Atom),
    /// `discriminant(x)`: decides, never stores.
    Discr(Atom),
    /// An index local: decides, never stores.
    Index(Local),
    /// Argument of the call whose result this is: decides unless the result
    /// only ever had its payload stored.
    CallArg(Atom),
    /// An argument of a call that operates on a `&mut` argument (insert,
    /// write, advance): the call failing is about the operation, not a
    /// verdict on any value handed to it. Never followed.
    CallArgMut,
    /// Argument of a call that could write through a `&mut` to the dest.
    ViaMut(Atom),
}

#[derive(Clone, Debug)]
struct Def {
    /// Field path of the destination place within its local.
    dest: Vec<u32>,
    reads: Vec<Read>,
}

/// A call to a crate-local fn: callee and argument atoms by position.
type LocalCall = (LocalDefId, Vec<Option<Atom>>);

struct Facts {
    defs: IndexVec<Local, Vec<Def>>,
    /// Whole-local moves (`_a = move _b`) and `Try::branch` look-through.
    alias: IndexVec<Local, Vec<Local>>,
    /// Locals whose variant payload was moved whole into this local.
    payload_alias: IndexVec<Local, Vec<Local>>,
    /// Atoms holding a `Self` value that reaches the return.
    self_roots: Vec<Atom>,
    local_calls: HashMap<Local, Vec<LocalCall>>,
    /// Blocks that assign a failure into the return, minus resource errors.
    failure_blocks: Vec<BasicBlock>,
    /// Closures defined in the body whose signature can yield `Self`.
    closures: Vec<LocalDefId>,
}

fn is_resource_error(tcx: TyCtxt<'_>, ty: Ty<'_>, extra: &[String]) -> bool {
    let ty::Adt(adt, _) = ty.peel_refs().kind() else {
        return false;
    };
    let entries = RESOURCE_ERRORS
        .iter()
        .copied()
        .chain(extra.iter().map(String::as_str));
    matches_config_path(tcx, adt.did(), entries)
}

/// Trait methods whose result is the receiver's value seen differently.
fn is_view_call(tcx: TyCtxt<'_>, callee: DefId) -> bool {
    let Some(tr) = tcx.trait_of_assoc(callee) else {
        return false;
    };
    let li = tcx.lang_items();
    [
        li.deref_trait(),
        li.deref_mut_trait(),
        li.index_trait(),
        li.index_mut_trait(),
        li.clone_trait(),
    ]
    .contains(&Some(tr))
        || [
            sym::Borrow,
            sym::BorrowMut,
            sym::AsRef,
            sym::AsMut,
            rustc_span::Symbol::intern("ToOwned"),
        ]
        .iter()
        .any(|s| tcx.is_diagnostic_item(*s, tr))
}

fn is_self_ty(ty: Ty<'_>, self_did: DefId) -> bool {
    matches!(ty.kind(), ty::Adt(a, _) if a.did() == self_did)
}

fn mentions_self(ty: Ty<'_>, self_did: DefId) -> bool {
    ty.walk()
        .any(|arg| arg.as_type().is_some_and(|t| is_self_ty(t, self_did)))
}

/// Reads of one operand/place in value position; `exact` tags the Exact case.
fn reads_of_place(place: Place<'_>, exact: impl FnOnce(Atom) -> Read, out: &mut Vec<Read>) {
    let info = place_info(place);
    out.extend(info.index_locals.iter().map(|l| Read::Index(*l)));
    out.push(match info.exactness {
        Exactness::VariantPayload => Read::Payload(info.atom),
        Exactness::Exact => exact(info.atom),
        Exactness::Inexact => Read::Derived(info.atom),
    });
}

fn reads_of_operand(op: &Operand<'_>, exact: impl FnOnce(Atom) -> Read, out: &mut Vec<Read>) {
    if let Some(p) = op.place() {
        reads_of_place(p, exact, out);
    }
}

/// `self_did` is the type under construction; None for a body whose failure
/// exits are the question and whose return value is not (no `Self` roots, no
/// closures to follow).
fn gather<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &Body<'tcx>,
    self_did: Option<DefId>,
    extra_resource_errors: &[String],
) -> Facts {
    let n = body.local_decls.len();
    let mut defs: IndexVec<Local, Vec<Def>> = IndexVec::from_elem_n(Vec::new(), n);
    let mut alias: IndexVec<Local, Vec<Local>> = IndexVec::from_elem_n(Vec::new(), n);
    let mut payload_alias: IndexVec<Local, Vec<Local>> = IndexVec::from_elem_n(Vec::new(), n);
    let mut mut_ref_to: HashMap<Local, Atom> = HashMap::new();
    let mut local_calls: HashMap<Local, Vec<LocalCall>> = HashMap::new();
    // Failure exits (Err/None with the error type) and `Self` payload roots,
    // each pending on whether their local reaches the return.
    let mut pending_failures: Vec<(Local, (BasicBlock, Option<Ty<'tcx>>))> = Vec::new();
    let mut pending_roots: Vec<(Local, Atom)> = Vec::new();
    let mut closures = Vec::new();

    for (bb, data) in body.basic_blocks.iter_enumerated() {
        if data.is_cleanup {
            continue;
        }
        for stmt in &data.statements {
            let StatementKind::Assign(assign) = &stmt.kind else {
                continue;
            };
            let (dest, rvalue) = &**assign;
            let dinfo = place_info(*dest);
            let mut reads = Vec::new();
            match rvalue {
                Rvalue::Use(op, _) => {
                    if dest.projection.is_empty()
                        && let Some(p) = op.place()
                    {
                        let pinfo = place_info(p);
                        if p.projection.is_empty() {
                            alias[dest.local].push(p.local);
                        } else if pinfo.exactness == Exactness::VariantPayload {
                            payload_alias[dest.local].push(p.local);
                        }
                    }
                    reads_of_operand(op, Read::Same, &mut reads);
                }
                Rvalue::Repeat(op, _)
                | Rvalue::Cast(_, op, _)
                | Rvalue::UnaryOp(_, op)
                | Rvalue::WrapUnsafeBinder(op, _) => {
                    reads_of_operand(op, Read::Derived, &mut reads)
                }
                Rvalue::Ref(_, _, place)
                | Rvalue::RawPtr(_, place)
                | Rvalue::Reborrow(_, _, place) => {
                    let is_mut = matches!(
                        rvalue,
                        Rvalue::Ref(_, BorrowKind::Mut { .. }, _)
                            | Rvalue::RawPtr(RawPtrKind::Mut, _)
                            | Rvalue::Reborrow(_, Mutability::Mut, _)
                    );
                    if is_mut && dest.projection.is_empty() {
                        mut_ref_to.insert(dest.local, place_info(*place).atom);
                    }
                    reads_of_place(*place, Read::Same, &mut reads);
                }
                Rvalue::CopyForDeref(place) => reads_of_place(*place, Read::Same, &mut reads),
                Rvalue::BinaryOp(_, ops) => {
                    let (a, b) = &**ops;
                    reads_of_operand(a, Read::Derived, &mut reads);
                    reads_of_operand(b, Read::Derived, &mut reads);
                }
                Rvalue::Discriminant(place) => {
                    let info = place_info(*place);
                    reads.extend(info.index_locals.iter().map(|l| Read::Index(*l)));
                    reads.push(Read::Discr(info.atom));
                }
                Rvalue::Aggregate(kind, ops) => {
                    let kind = &**kind;
                    for (i, op) in ops.iter().enumerate() {
                        reads_of_operand(op, |a| Read::AggField(i as u32, a), &mut reads);
                    }
                    match kind {
                        AggregateKind::Adt(did, vidx, _, _, None)
                            if dest.projection.is_empty()
                                && (tcx.is_diagnostic_item(sym::Result, *did)
                                    || tcx.is_diagnostic_item(sym::Option, *did)) =>
                        {
                            let first = ops.iter().next();
                            let name = tcx.adt_def(*did).variant(*vidx).name;
                            if name == sym::Err || name == sym::None {
                                let err = first.map(|o| o.ty(&body.local_decls, tcx));
                                pending_failures.push((dest.local, (bb, err)));
                            } else if let Some(p) = first.and_then(Operand::place) {
                                pending_roots.push((dest.local, place_info(p).atom));
                            }
                        }
                        AggregateKind::Closure(cdid, args) => {
                            let out = args.as_closure().sig().output().skip_binder();
                            if self_did.is_some_and(|s| mentions_self(out, s))
                                && let Some(l) = cdid.as_local()
                            {
                                closures.push(l);
                            }
                        }
                        _ => {}
                    }
                }
                Rvalue::ThreadLocalRef(_) => {}
            }
            defs[dest.local].push(Def {
                dest: dinfo.atom.path,
                reads,
            });
        }

        let Some(term) = &data.terminator else {
            continue;
        };
        if let TerminatorKind::Call {
            func,
            args,
            destination,
            ..
        } = &term.kind
        {
            let arg_atoms: Vec<Option<Atom>> = args
                .iter()
                .map(|a| a.node.place().map(|p| place_info(p).atom))
                .collect();
            // A call handed a `&mut` is an operation on that argument (insert,
            // write, advance); its other arguments are the data operated
            // with, and its failing judges none of them.
            let operation = arg_atoms
                .iter()
                .flatten()
                .any(|a| a.path.is_empty() && mut_ref_to.contains_key(&a.local));
            let call_arg = |a: &Atom| {
                if operation {
                    Read::CallArgMut
                } else {
                    Read::CallArg(a.clone())
                }
            };
            let mut reads: Vec<Read> = arg_atoms.iter().flatten().map(call_arg).collect();
            // A `&mut x` handed to the call: whatever the callee decides from
            // its arguments can land in `x`.
            for a in arg_atoms.iter().flatten() {
                if a.path.is_empty()
                    && let Some(target) = mut_ref_to.get(&a.local)
                {
                    defs[target.local].push(Def {
                        dest: target.path.clone(),
                        reads: arg_atoms
                            .iter()
                            .flatten()
                            .cloned()
                            .map(Read::ViaMut)
                            .collect(),
                    });
                }
            }
            if let Some((callee, cargs)) = func.const_fn_def() {
                if is_view_call(tcx, callee)
                    && let Some(Some(recv)) = arg_atoms.first()
                {
                    // `deref`, `index`, `clone`, `as_ref`, `borrow`: the result
                    // is (a view of) the receiver's value, not a verdict on it.
                    reads.clear();
                    let li = tcx.lang_items();
                    let tr = tcx.trait_of_assoc(callee);
                    if tr.is_some() && (tr == li.index_trait() || tr == li.index_mut_trait()) {
                        reads.push(Read::Same(recv.extended(&[ANY_ELEM])));
                    } else {
                        reads.push(Read::Same(recv.clone()));
                    }
                    for a in arg_atoms.iter().skip(1).flatten() {
                        if a.path.is_empty() {
                            reads.push(Read::Index(a.local));
                        }
                    }
                } else if tcx.is_lang_item(callee, LangItem::TryTraitFromResidual) {
                    let residual = cargs
                        .get(1)
                        .and_then(|g| g.as_type())
                        .or_else(|| args.first().map(|a| a.node.ty(&body.local_decls, tcx)));
                    let err = residual.and_then(|r| result_err_ty(tcx, r));
                    pending_failures.push((destination.local, (bb, err)));
                } else if tcx.is_lang_item(callee, LangItem::TryTraitBranch) {
                    if let Some(Some(r)) = arg_atoms.first()
                        && r.path.is_empty()
                    {
                        alias[destination.local].push(r.local);
                    }
                } else if let Some(l) = callee.as_local()
                    && destination.projection.is_empty()
                {
                    local_calls
                        .entry(destination.local)
                        .or_default()
                        .push((l, arg_atoms.clone()));
                    reads.clear();
                    reads.extend(arg_atoms.iter().flatten().map(call_arg));
                }
            }
            let dinfo = place_info(*destination);
            defs[destination.local].push(Def {
                dest: dinfo.atom.path,
                reads,
            });
        }
    }

    let returned = closure_over(RETURN_PLACE, |l| alias[l].clone());

    let mut failure_blocks: Vec<BasicBlock> = reaching(pending_failures, &returned)
        .into_iter()
        .filter(|(_, err)| !err.is_some_and(|t| is_resource_error(tcx, t, extra_resource_errors)))
        .map(|(bb, _)| bb)
        .collect();
    let mut self_roots = reaching(pending_roots, &returned);
    // A body returning bare `Self` (a helper, a closure): the return place
    // itself is the self value.
    if self_did.is_some_and(|s| is_self_ty(body.local_decls[RETURN_PLACE].ty, s)) {
        self_roots.push(Atom::whole(RETURN_PLACE));
    }
    failure_blocks.sort();
    failure_blocks.dedup();

    Facts {
        defs,
        alias,
        payload_alias,
        self_roots,
        local_calls,
        failure_blocks,
        closures,
    }
}

/// The items whose local reaches the return.
fn reaching<T>(pending: Vec<(Local, T)>, returned: &HashSet<Local>) -> Vec<T> {
    let kept = pending.into_iter().filter(|(l, _)| returned.contains(l));
    kept.map(|(_, t)| t).collect()
}

fn closure_over(start: Local, mut succ: impl FnMut(Local) -> Vec<Local>) -> HashSet<Local> {
    let mut seen = HashSet::new();
    let mut q = VecDeque::from([start]);
    while let Some(l) = q.pop_front() {
        if seen.insert(l) {
            q.extend(succ(l));
        }
    }
    seen
}

// ── slices ───────────────────────────────────────────────────────────────────

/// Storage slice from `roots`: every atom whose value ends up (as itself, or
/// arithmetically combined) in a root. Returns the slice and the atoms whose
/// variant payload was consumed on the way (produced, not inspected).
fn storage_slice(defs: &IndexVec<Local, Vec<Def>>, roots: &[Atom]) -> (Vec<Atom>, Vec<Atom>) {
    let mut produced: Vec<Atom> = Vec::new();
    let slice = walk_slice(
        defs,
        roots,
        |_| (),
        |(), r| match r {
            Read::Derived(a) => Some(a.clone()),
            Read::Payload(a) => {
                produced.push(a.clone());
                None
            }
            _ => None,
        },
    );
    (slice, produced)
}

/// Decision slice from a branch operand. `opaque` holds atoms whose defining
/// call must not be entered: produced-payload bases (the branch consumed the
/// callee's verdict, it did not inspect the arguments) and stored values (a
/// method failing on a stored value implicates that value, not everything
/// its constructor was handed).
fn decision_slice(defs: &IndexVec<Local, Vec<Def>>, roots: &[Atom], opaque: &[Atom]) -> Vec<Atom> {
    walk_slice(
        defs,
        roots,
        |at| !opaque.iter().any(|p| p.overlaps(at)),
        |&enter_calls, r| match r {
            Read::Derived(a) | Read::Payload(a) | Read::Discr(a) | Read::ViaMut(a) => {
                Some(a.clone())
            }
            Read::Index(l) => Some(Atom::whole(*l)),
            Read::CallArg(a) => enter_calls.then(|| a.clone()),
            Read::Same(_) | Read::AggField(..) | Read::CallArgMut => None,
        },
    )
}

/// Worklist closure of `roots` over `defs`. `Same` and `AggField` compose
/// paths here; any other read is followed to whatever atom `other` names,
/// given the per-atom state `enter` computed once when the atom was dequeued.
fn walk_slice<S>(
    defs: &IndexVec<Local, Vec<Def>>,
    roots: &[Atom],
    mut enter: impl FnMut(&Atom) -> S,
    mut other: impl FnMut(&S, &Read) -> Option<Atom>,
) -> Vec<Atom> {
    let mut seen: HashSet<Atom> = HashSet::new();
    let mut q: VecDeque<Atom> = roots.iter().cloned().collect();
    while let Some(at) = q.pop_front() {
        if !seen.insert(at.clone()) {
            continue;
        }
        let s = enter(&at);
        for def in &defs[at.local] {
            let (below, rem): (bool, &[u32]) = if at.path.starts_with(&def.dest) {
                (true, &at.path[def.dest.len()..])
            } else if def.dest.starts_with(&at.path) {
                (false, &[])
            } else {
                continue;
            };
            for r in &def.reads {
                q.extend(match r {
                    Read::Same(a) => Some(if below { a.extended(rem) } else { a.clone() }),
                    Read::AggField(i, a) if below && !rem.is_empty() => {
                        (rem[0] == *i).then(|| a.extended(&rem[1..]))
                    }
                    Read::AggField(_, a) => Some(a.clone()),
                    Read::Derived(_)
                    | Read::Payload(_)
                    | Read::Discr(_)
                    | Read::Index(_)
                    | Read::CallArg(_)
                    | Read::CallArgMut
                    | Read::ViaMut(_) => other(&s, r),
                });
            }
        }
    }
    seen.into_iter().collect()
}

fn slices_meet(storage: &[Atom], decision: &[Atom]) -> bool {
    let mut by_local: HashMap<Local, Vec<&Atom>> = HashMap::new();
    for d in decision {
        by_local.entry(d.local).or_default().push(d);
    }
    storage.iter().any(|s| {
        by_local
            .get(&s.local)
            .is_some_and(|ds| ds.iter().any(|d| d.inspects(s)))
    })
}

// ── control dependence ───────────────────────────────────────────────────────

/// The branch at `bb` switches on the discriminant of a `Result<_, E>` whose
/// `E` is a resource error: its outcome is the environment's, whatever the
/// constructor then returns.
fn decides_on_resource<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &Body<'tcx>,
    facts: &Facts,
    bb: BasicBlock,
    extra: &[String],
) -> bool {
    let TerminatorKind::SwitchInt { discr, .. } = &body.basic_blocks[bb].terminator().kind else {
        return false;
    };
    let Some(p) = discr.place() else {
        return false;
    };
    // `_d = discriminant(_r); switchInt(move _d)`: look one def back.
    let mut candidates = vec![p.local];
    for def in &facts.defs[p.local] {
        for r in &def.reads {
            if let Read::Discr(a) = r {
                candidates.push(a.local);
                // `?`: `_t = Try::branch(_r)`; the residual type sits on `_r`.
                candidates.extend(facts.alias[a.local].iter().copied());
            }
        }
    }
    candidates.into_iter().any(|l| {
        result_err_ty(tcx, body.local_decls[l].ty).is_some_and(|e| is_resource_error(tcx, e, extra))
    })
}

// ── per-field roots, helpers, and the verdict ────────────────────────────────

type Summary = HashMap<FieldIdx, Vec<usize>>;
type Memo = HashMap<LocalDefId, Option<Summary>>;

/// Which parameters of helper `def` (returning `Self`, `Result<Self, _>` or
/// `Option<Self>`) flow into which field of the `Self` it builds.
fn helper_summary<'tcx>(
    tcx: TyCtxt<'tcx>,
    def: LocalDefId,
    self_did: DefId,
    memo: &mut Memo,
    depth: usize,
) -> Option<Summary> {
    if let Some(m) = memo.get(&def) {
        return m.clone();
    }
    if depth > 3 {
        return None;
    }
    memo.insert(def, None);
    let body = mir_for(tcx, def)?;
    let facts = gather(tcx, &body, Some(self_did), &[]);
    let nfields = tcx.adt_def(self_did).non_enum_variant().fields.len();
    let roots = field_roots(tcx, &body, &facts, self_did, nfields, memo, depth + 1);
    let argc = body.arg_count;
    let mut out: Summary = HashMap::new();
    for (f, r) in roots {
        let (slice, _) = storage_slice(&facts.defs, &r);
        let params: Vec<usize> = (1..=argc)
            .filter(|i| slice.iter().any(|a| a.local.as_usize() == *i))
            .map(|i| i - 1)
            .collect();
        if !params.is_empty() {
            out.entry(f).or_default().extend(params);
        }
    }
    memo.insert(def, Some(out.clone()));
    Some(out)
}

fn root_fields(roots: &mut HashMap<FieldIdx, Vec<Atom>>, base: &Atom, nfields: usize) {
    for f in 0..nfields {
        roots
            .entry(FieldIdx::from_usize(f))
            .or_default()
            .push(base.extended(&[f as u32]));
    }
}

/// field -> root atoms: `carrier.f` for every `Self`-typed local the returned
/// value passes through, plus helper-call arguments the helper stores in `f`.
fn field_roots<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &Body<'tcx>,
    facts: &Facts,
    self_did: DefId,
    nfields: usize,
    memo: &mut Memo,
    depth: usize,
) -> HashMap<FieldIdx, Vec<Atom>> {
    let mut roots: HashMap<FieldIdx, Vec<Atom>> = HashMap::new();
    let mut carriers: HashSet<Local> = HashSet::new();
    for r in &facts.self_roots {
        if r.path.is_empty() {
            carriers.extend(closure_over(r.local, |l| {
                let mut v = facts.alias[l].clone();
                v.extend(facts.payload_alias[l].iter().copied());
                v
            }));
        } else {
            // `Ok(pair.1)`-style: the self value is a sub-place; root there.
            root_fields(&mut roots, r, nfields);
        }
    }
    for &c in &carriers {
        if is_self_ty(body.local_decls[c].ty, self_did) {
            root_fields(&mut roots, &Atom::whole(c), nfields);
        }
        if let Some(calls) = facts.local_calls.get(&c) {
            for (callee, args) in calls {
                let ret = tcx
                    .fn_sig(callee.to_def_id())
                    .instantiate_identity()
                    .skip_normalization()
                    .output()
                    .skip_binder();
                if !mentions_self(ret, self_did) {
                    continue;
                }
                if let Some(summary) = helper_summary(tcx, *callee, self_did, memo, depth) {
                    for (f, params) in summary {
                        for p in params {
                            if let Some(Some(a)) = args.get(p) {
                                roots.entry(f).or_default().push(a.clone());
                            }
                        }
                    }
                }
            }
        }
    }
    roots
}

fn analyze_body<'tcx>(
    tcx: TyCtxt<'tcx>,
    def: LocalDefId,
    self_did: DefId,
    extra_resource_errors: &[String],
    memo: &mut Memo,
    depth: usize,
    out: &mut HashMap<FieldIdx, Span>,
) {
    if depth > 3 {
        return;
    }
    let Some(body) = mir_for(tcx, def) else {
        return;
    };
    if body.tainted_by_errors.is_some() {
        return;
    }
    if let Some(e) = result_err_ty(tcx, body.local_decls[RETURN_PLACE].ty)
        && is_resource_error(tcx, e, extra_resource_errors)
    {
        return;
    }
    let facts = gather(tcx, &body, Some(self_did), extra_resource_errors);
    let closures = facts.closures.clone();
    if !facts.failure_blocks.is_empty() {
        let nfields = tcx.adt_def(self_did).non_enum_variant().fields.len();
        let roots = field_roots(tcx, &body, &facts, self_did, nfields, memo, depth);
        let mut slices: Vec<(FieldIdx, Vec<Atom>)> = Vec::new();
        let mut opaque: Vec<Atom> = Vec::new();
        for (f, r) in roots {
            let (s, p) = storage_slice(&facts.defs, &r);
            opaque.extend(p);
            opaque.extend(s.iter().cloned());
            slices.push((f, s));
        }
        if !slices.is_empty() {
            let flow = FlowGraph::new(&body);
            let mut switch_slices: HashMap<BasicBlock, Vec<Atom>> = HashMap::new();
            for &fb in &facts.failure_blocks {
                for branch in flow.control_deps(fb) {
                    if decides_on_resource(tcx, &body, &facts, branch, extra_resource_errors) {
                        continue;
                    }
                    let ds = switch_slices.entry(branch).or_insert_with(|| {
                        decision_slice(&facts.defs, &switch_operand_atoms(&body, branch), &opaque)
                    });
                    for (f, ss) in &slices {
                        if !out.contains_key(f) && slices_meet(ss, ds) {
                            out.insert(*f, body.basic_blocks[branch].terminator().source_info.span);
                        }
                    }
                }
            }
        }
    }
    drop(body);
    // A closure that can yield `Self` (`.and_then(|v| ...)`, `try_parse(|i|
    // ...)`) is a constructor body of its own.
    for c in closures {
        analyze_body(
            tcx,
            c,
            self_did,
            extra_resource_errors,
            memo,
            depth + 1,
            out,
        );
    }
}
