//! Which variants of an enum can a function return? MIR-level and
//! conservative: every value reaching the return place is traced back through
//! plain copies between locals to either an enum aggregate with a known
//! variant, a constant whose variant the valtree names, or something
//! untraceable — a parameter, a call result, a projection — in which case the
//! whole set is unknowable and the answer is `None`.
//!
//! The tracing is flow-insensitive per local: a local's possible sources are
//! everything ever assigned to it. That over-approximates the returned set
//! (an overwritten variant still counts), which for the narrowing lints is
//! the safe direction — a variant wrongly *in* the set can only silence a
//! finding, never invent one.

use std::collections::{HashMap, HashSet};

use rustc_abi::VariantIdx;
use rustc_hir::def_id::DefId;
use rustc_lint::LateContext;
use rustc_middle::mir::{
    AggregateKind, Local, Operand, RETURN_PLACE, Rvalue, StatementKind, TerminatorKind,
};
use rustc_middle::ty;
use rustc_span::def_id::LocalDefId;

use crate::mir_flow::mir_for;

enum Source {
    Variant(VariantIdx),
    Alias(Local),
    Unknown,
}

/// The complete set of variants of `enum_did` that `func` can return, or
/// `None` when any value reaching the return is untraceable.
pub(crate) fn returned_variants(
    cx: &LateContext<'_>,
    func: LocalDefId,
    enum_did: DefId,
) -> Option<HashSet<VariantIdx>> {
    let body = mir_for(cx.tcx, func)?;
    let mut sources: HashMap<Local, Vec<Source>> = HashMap::new();

    for block in body.basic_blocks.iter() {
        for stmt in &block.statements {
            let StatementKind::Assign(assign) = &stmt.kind else {
                continue;
            };
            let (place, rvalue) = &**assign;
            // Projections on the destination (field stores) never build the
            // whole enum value; ignore them rather than misattribute.
            if !place.projection.is_empty() {
                continue;
            }
            let src = match rvalue {
                Rvalue::Aggregate(kind, _) => match &**kind {
                    AggregateKind::Adt(did, vidx, ..) if *did == enum_did => Source::Variant(*vidx),
                    _ => Source::Unknown,
                },
                Rvalue::Use(op, _) => operand_source(cx, op, enum_did),
                _ => Source::Unknown,
            };
            sources.entry(place.local).or_default().push(src);
        }
        // A call writing its result into a local makes that local
        // untraceable; so does anything else a terminator can write.
        if let TerminatorKind::Call { destination, .. } = &block.terminator().kind {
            sources
                .entry(destination.local)
                .or_default()
                .push(Source::Unknown);
        }
    }

    // Resolve the return place through aliases.
    let mut set = HashSet::new();
    let mut stack = vec![RETURN_PLACE];
    let mut seen: HashSet<Local> = HashSet::new();
    while let Some(local) = stack.pop() {
        if !seen.insert(local) {
            continue;
        }
        // A local with no recorded assignment is an argument or untouched
        // storage: untraceable.
        let srcs = sources.get(&local)?;
        for src in srcs {
            match src {
                Source::Variant(v) => {
                    set.insert(*v);
                }
                Source::Alias(l) => stack.push(*l),
                Source::Unknown => return None,
            }
        }
    }
    (!set.is_empty()).then_some(set)
}

/// What an operand contributes: a plain local copy aliases it; a constant of
/// the enum type names its variant through the valtree; anything else is
/// untraceable.
fn operand_source<'tcx>(cx: &LateContext<'tcx>, op: &Operand<'tcx>, enum_did: DefId) -> Source {
    match op {
        Operand::Copy(p) | Operand::Move(p) => {
            if p.projection.is_empty() {
                Source::Alias(p.local)
            } else {
                Source::Unknown
            }
        }
        // Whatever future operand kinds carry, they are not traceable
        // enum values built here.
        Operand::RuntimeChecks(_) => Source::Unknown,
        Operand::Constant(c) => {
            let ty = c.const_.ty();
            let ty::Adt(adt, _) = ty.kind() else {
                return Source::Unknown;
            };
            if adt.did() != enum_did {
                return Source::Unknown;
            }
            // Best effort: when the constant evaluates to a scalar equal to
            // exactly one variant's discriminant, that variant is the source.
            // Anything else (niche encodings, by-ref payloads) is unknowable.
            let Some(scalar) = c.const_.try_eval_scalar_int(cx.tcx, cx.typing_env()) else {
                return Source::Unknown;
            };
            let val = scalar.to_uint(scalar.size());
            let mut hit: Option<VariantIdx> = None;
            for (vidx, discr) in adt.discriminants(cx.tcx) {
                if discr.val == val {
                    if hit.is_some() {
                        return Source::Unknown;
                    }
                    hit = Some(vidx);
                }
            }
            match hit {
                Some(v) => Source::Variant(v),
                None => Source::Unknown,
            }
        }
    }
}
