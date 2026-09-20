//! Builds the text of the finding: how each value is named, the message, the
//! notes and the help lines.

use clippy_utils::source::snippet_opt;
use rustc_data_structures::fx::FxHashMap;
use rustc_hir::attrs::InlineAttr;
use rustc_hir::def_id::LocalDefId;
use rustc_index::IndexVec;
use rustc_index::bit_set::DenseBitSet;
use rustc_lint::LateContext;
use rustc_middle::mir::visit::Visitor;
use rustc_middle::mir::{
    self, BasicBlock, Local, Location, Operand, Place, RETURN_PLACE, Rvalue, Statement,
    StatementKind, TerminatorKind,
};
use rustc_middle::ty::print::with_no_trimmed_paths;
use rustc_middle::ty::{self, GenericArgsRef, GenericParamDefKind, Ty, TyCtxt, TypeVisitableExt};
use rustc_span::{Span, Symbol};

use super::GENERIC_BODY_NOT_GENERIC;
use super::classify::{PlaceParams, counted_statement, counted_terminator};
use super::source_span::SourceSpan;
use crate::baseline::{emit_hir_then, join};
use crate::mir_flow::{local_name, reads_any};

fn inline_note(tcx: TyCtxt<'_>, def: LocalDefId) -> Option<String> {
    matches!(tcx.codegen_fn_attrs(def).inline, InlineAttr::Hint).then(|| {
        format!(
            "`{}` is `#[inline]`. Leave `#[inline]` off the separate function, or every \
             codegen unit that calls it compiles its own copy and shares nothing",
            tcx.def_path_str(def)
        )
    })
}

/// `bytes.as_ref()` giving a `&[u8]`: `arg` through `AsRef::as_ref` obtains `&[u8]`.
pub(super) struct Conversion<'tcx> {
    arg: mir::Local,
    through: String,
    obtains: Ty<'tcx>,
}

/// The blocks every call runs exactly once, in order.
fn entry_blocks(body: &mir::Body<'_>) -> Vec<BasicBlock> {
    let predecessors = body.basic_blocks.predecessors();
    let mut run = Vec::new();
    if !predecessors[mir::START_BLOCK].is_empty() {
        return run;
    }
    let mut block = mir::START_BLOCK;
    while run.len() <= body.basic_blocks.len() {
        run.push(block);
        let next = match body.basic_blocks[block].terminator().kind {
            TerminatorKind::Goto { target }
            | TerminatorKind::Drop { target, .. }
            | TerminatorKind::Assert { target, .. }
            | TerminatorKind::Call {
                target: Some(target),
                ..
            } => target,
            _ => break,
        };
        if predecessors[next].len() != 1 || body.basic_blocks[next].is_cleanup {
            break;
        }
        block = next;
    }
    run
}

/// Non-empty when every dependent item is a trait-method call on a generic argument in
/// `entry_blocks`, a borrow or move into one, or a drop. The caller checks `owns_signature`.
pub(super) fn only_conversions<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
) -> Vec<Conversion<'tcx>> {
    let decls = &body.local_decls;
    let bare_param_arg = |local: mir::Local| {
        body.local_kind(local) == mir::LocalKind::Arg
            && matches!(decls[local].ty.peel_refs().kind(), ty::Param(_))
    };
    if !body.args_iter().any(bare_param_arg) {
        return Vec::new();
    }
    let mut arg_of_temp: FxHashMap<mir::Local, mir::Local> = FxHashMap::default();
    // Such an argument or a temporary holding one, as is or dereferenced.
    let arg_of_place =
        |temps: &FxHashMap<mir::Local, mir::Local>, place: Place<'tcx>| -> Option<mir::Local> {
            if !place
                .projection
                .iter()
                .all(|elem| matches!(elem, mir::ProjectionElem::Deref))
            {
                return None;
            }
            if bare_param_arg(place.local) {
                Some(place.local)
            } else {
                temps.get(&place.local).copied()
            }
        };
    let run = entry_blocks(body);
    let mut in_entry_blocks: IndexVec<BasicBlock, bool> =
        IndexVec::from_elem_n(false, body.basic_blocks.len());
    for &block in &run {
        in_entry_blocks[block] = true;
    }
    let rest = body
        .basic_blocks
        .indices()
        .filter(|&block| !in_entry_blocks[block] && !body.basic_blocks[block].is_cleanup);
    let order: Vec<BasicBlock> = run.iter().copied().chain(rest).collect();
    let mut places = PlaceParams::new(tcx, body);
    let mut found: Vec<Conversion<'tcx>> = Vec::new();
    for block in order {
        let data = &body.basic_blocks[block];
        let early = in_entry_blocks[block];
        for statement in &data.statements {
            if !counted_statement(statement) || !places.dependent_statement(statement) {
                continue;
            }
            if early
                && let Some((target, rvalue)) = statement.kind.as_assign()
                && let Some(temp) = target.as_local()
                && body.local_kind(temp) == mir::LocalKind::Temp
                && let mir::Rvalue::Ref(_, mir::BorrowKind::Shared, place)
                | mir::Rvalue::CopyForDeref(place)
                | mir::Rvalue::Use(mir::Operand::Move(place) | mir::Operand::Copy(place), _) =
                    rvalue
                && let Some(arg) = arg_of_place(&arg_of_temp, *place)
            {
                arg_of_temp.insert(temp, arg);
            } else {
                return Vec::new();
            }
        }
        let Some(terminator) = &data.terminator else {
            continue;
        };
        if !counted_terminator(terminator) || !places.dependent_terminator(terminator) {
            continue;
        }
        match &terminator.kind {
            TerminatorKind::Drop { place, .. }
                if place
                    .as_local()
                    .is_some_and(|l| bare_param_arg(l) || arg_of_temp.contains_key(&l)) => {}
            TerminatorKind::Call {
                func,
                args,
                destination,
                target: Some(_),
                ..
            } if early => {
                if let [operand] = &args[..]
                    && let Some(place) = operand.node.place()
                    && let Some(arg) = arg_of_place(&arg_of_temp, place)
                    && let Some((callee, _)) = func.const_fn_def()
                    && let Some(trait_id) = tcx.trait_of_assoc(callee)
                    && let Some(obtained) = destination.as_local()
                    && !decls[obtained].ty.has_non_region_param()
                    && !decls[obtained].ty.is_unit()
                    && !decls[obtained].ty.is_never()
                {
                    found.push(Conversion {
                        arg,
                        // `item_name` is safe before a diagnostic is certain.
                        through: format!("{}::{}", tcx.item_name(trait_id), tcx.item_name(callee)),
                        obtains: decls[obtained].ty,
                    });
                } else {
                    return Vec::new();
                }
            }
            _ => return Vec::new(),
        }
    }
    if body
        .args_iter()
        .filter(|&arg| bare_param_arg(arg))
        .any(|arg| !found.iter().any(|c| c.arg == arg))
    {
        return Vec::new();
    }
    found
}

/// `` `bytes` ``, or "its `<type>` argument" for an unnamed pattern argument.
fn arg_name(body: &mir::Body<'_>, arg: mir::Local) -> String {
    match local_name(body, arg) {
        Some(name) => format!("`{name}`"),
        None => with_no_trimmed_paths!(format!("its `{}` argument", body.local_decls[arg].ty)),
    }
}

/// The second help line without the fn's name: `def_path_str` may only run
/// once a diagnostic is certain, so `report` prepends it.
pub(super) fn conversion_help(
    body: &mir::Body<'_>,
    conversions: &[Conversion<'_>],
) -> Option<String> {
    let mut clauses: Vec<String> = Vec::new();
    let mut listed = 0;
    let mut single: Option<String> = None;
    for arg in body.args_iter() {
        let mut ways: Vec<String> = Vec::new();
        for conversion in conversions.iter().filter(|c| c.arg == arg) {
            let ty = with_no_trimmed_paths!(format!("`{}`", conversion.obtains));
            let way = format!("a {ty} through `{}`", conversion.through);
            if !ways.contains(&way) {
                ways.push(way);
                listed += 1;
                single = Some(ty);
            }
        }
        if !ways.is_empty() {
            clauses.push(format!(
                "{} only to get {}",
                arg_name(body, arg),
                join(&ways, "and")
            ));
        }
    }
    let instead = match (listed, single) {
        (1, Some(one)) => one,
        (0, _) | (_, None) => return None,
        _ => "those".to_owned(),
    };
    Some(format!(
        "uses {}, so it could take {instead} instead",
        join(&clauses, "and")
    ))
}

fn param_names(tcx: TyCtxt<'_>, def: LocalDefId) -> Vec<String> {
    let generics = tcx.generics_of(def);
    (0..generics.count())
        .map(|i| generics.param_at(i, tcx))
        .filter(|p| !matches!(p.kind, GenericParamDefKind::Lifetime))
        .map(|p| format!("`{}`", p.name))
        .collect()
}

fn render_args<'tcx>(tcx: TyCtxt<'tcx>, def: LocalDefId, args: GenericArgsRef<'tcx>) -> String {
    let generics = tcx.generics_of(def);
    let pairs: Vec<String> = args
        .iter()
        .enumerate()
        .filter_map(|(i, arg)| {
            let param = generics.param_at(i, tcx);
            (!matches!(param.kind, GenericParamDefKind::Lifetime))
                .then(|| with_no_trimmed_paths!(format!("{} = {arg}", param.name)))
        })
        .collect();
    pairs.join(", ")
}

enum AssignedTo {
    Let(Symbol),
    Field(Symbol),
}

fn field_name<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    place: Place<'tcx>,
) -> Option<Symbol> {
    let (base, mir::ProjectionElem::Field(field, _)) = place.as_ref().last_projection()? else {
        return None;
    };
    let base_ty = base.ty(body, tcx);
    let ty::Adt(def, _) = *base_ty.ty.kind() else {
        return None;
    };
    let variant = match base_ty.variant_index {
        Some(index) => def.variant(index),
        None if def.is_enum() => return None,
        None => def.non_enum_variant(),
    };
    Some(variant.fields.get(field)?.name)
}

/// Where the one statement in `among` reading `local` assigns it, whole.
fn assigned_to<'tcx>(
    tcx: TyCtxt<'tcx>,
    body: &mir::Body<'tcx>,
    among: &DenseBitSet<BasicBlock>,
    local: Local,
) -> Option<AssignedTo> {
    let mut of = DenseBitSet::new_empty(body.local_decls.len());
    of.insert(local);
    let mut user: Option<&Statement<'tcx>> = None;
    for block in among.iter() {
        let data = &body.basic_blocks[block];
        if data.is_cleanup {
            continue;
        }
        for statement in &data.statements {
            let found = reads_any(&of, |uses| uses.visit_statement(statement, Location::START));
            let writes = match &statement.kind {
                StatementKind::Assign(assign) => assign.0.as_local() == Some(local),
                _ => false,
            };
            if found && !writes {
                if user.is_some() {
                    return None;
                }
                user = Some(statement);
            }
        }
        if let Some(terminator) = &data.terminator {
            let found = reads_any(&of, |uses| {
                uses.visit_terminator(terminator, Location::START)
            });
            let writes = matches!(
                &terminator.kind,
                TerminatorKind::Call { destination, .. } if destination.local == local
            );
            if found && !writes {
                return None;
            }
        }
    }
    let StatementKind::Assign(assign) = &user?.kind else {
        return None;
    };
    let (dest, rvalue) = &**assign;
    let whole = |operand: &Operand<'tcx>| match operand {
        Operand::Copy(place) | Operand::Move(place) => place.as_local() == Some(local),
        _ => false,
    };
    match rvalue {
        Rvalue::Use(operand, _) if whole(operand) => match dest.as_local() {
            Some(var) => local_name(body, var).map(AssignedTo::Let),
            None => field_name(tcx, body, *dest).map(AssignedTo::Field),
        },
        Rvalue::Ref(_, _, place)
            if place.local == local
                && matches!(place.projection[..], [mir::ProjectionElem::Deref]) =>
        {
            let var = dest.as_local()?;
            (body.local_decls[var].ty == body.local_decls[local].ty)
                .then(|| local_name(body, var).map(AssignedTo::Let))?
        }
        Rvalue::Aggregate(kind, operands) => {
            let mir::AggregateKind::Adt(adt, variant, _, _, None) = **kind else {
                return None;
            };
            let index = operands.iter().position(whole)?;
            let field = tcx
                .adt_def(adt)
                .variant(variant)
                .fields
                .get(rustc_abi::FieldIdx::from_usize(index))?;
            Some(AssignedTo::Field(field.name))
        }
        _ => None,
    }
}

/// `` `acc: u32` ``, "the returned `u32`", "the `u32` assigned to `total`",
/// "the `usize` from `src.len()`". Untrimmed: trimming needs a diagnostic.
fn render_local<'tcx>(
    cx: &LateContext<'tcx>,
    body: &mir::Body<'tcx>,
    body_span: Span,
    among: &DenseBitSet<BasicBlock>,
    local: Local,
) -> String {
    let decl = &body.local_decls[local];
    let ty = with_no_trimmed_paths!(decl.ty.to_string());
    if let Some(name) = local_name(body, local) {
        return format!("`{name}: {ty}`");
    }
    if local == RETURN_PLACE {
        return format!("the returned `{ty}`");
    }
    match assigned_to(cx.tcx, body, among, local) {
        Some(AssignedTo::Let(name)) => return format!("the `{ty}` assigned to `{name}`"),
        Some(AssignedTo::Field(name)) => return format!("the `{ty}` for field `{name}`"),
        None => {}
    }
    // A temporary's span is its expression. An argument's is outside the body.
    let span = decl.source_info.span;
    if !span.from_expansion()
        && body_span.contains(span)
        && let Some(text) = snippet_opt(cx, span)
        && !text.contains('\n')
        && text.len() <= 40
    {
        return format!("the `{ty}` from `{text}`");
    }
    format!("a `{ty}` temporary")
}

/// Skips `()` (MIR's `return` reads `_0` even in a unit fn) and repeats.
pub(super) fn render_locals<'tcx>(
    cx: &LateContext<'tcx>,
    def: LocalDefId,
    body: &mir::Body<'tcx>,
    locals: &[Local],
    among: &DenseBitSet<BasicBlock>,
) -> Vec<String> {
    let body_span = cx.tcx.hir_body_owned_by(def).value.span;
    let mut out: Vec<String> = Vec::new();
    for &local in locals {
        if body.local_decls[local].ty.is_unit() {
            continue;
        }
        let shown = render_local(cx, body, body_span, among, local);
        if !out.contains(&shown) {
            out.push(shown);
        }
    }
    out
}

fn finding_message(
    tcx: TyCtxt<'_>,
    def: LocalDefId,
    sets: usize,
    size: usize,
    total: usize,
) -> String {
    let name = tcx.def_path_str(def);
    let params = param_names(tcx, def);
    format!(
        "`{name}` is generic over {} and is compiled {sets} times in this crate, but one part of \
         it, {size} of its {total} statements (as MIR), does not use {}",
        join(&params, "and"),
        join(&params, "or"),
    )
}

fn part_note(part: &Finding, min_statements: usize) -> String {
    let reads = match &part.reads[..] {
        [] => "read nothing computed before them".to_owned(),
        list => format!("read only {}", list.join(", ")),
    };
    let produces = match &part.produces[..] {
        [] => "produce nothing read after them".to_owned(),
        list => format!("produce {}", list.join(", ")),
    };
    let mut rest = format!("{reads} and {produces}");
    match part.other_parts {
        0 => {}
        1 => rest.push_str(&format!(
            ", and the body has 1 more such part of at least {min_statements} statements"
        )),
        more => rest.push_str(&format!(
            ", and the body has {more} more such parts of at least {min_statements} statements"
        )),
    }
    rest
}

/// Rendered strings only, so the MIR body's borrow ends early. `total` is
/// the body's counted items.
pub(super) struct Finding {
    pub(super) def: LocalDefId,
    pub(super) site: Option<SourceSpan>,
    pub(super) size: usize,
    pub(super) reads: Vec<String>,
    pub(super) produces: Vec<String>,
    pub(super) other_parts: usize,
    pub(super) total: usize,
    pub(super) signature_help: Option<String>,
}

/// Emits with the fn's HirId so an `#[allow]` on the fn is honoured: during
/// `check_crate_post` the current node is the crate root.
pub(super) fn report<'tcx>(
    cx: &LateContext<'tcx>,
    part: &Finding,
    sets: usize,
    site: Option<(Span, GenericArgsRef<'tcx>)>,
    min_statements: usize,
) {
    let tcx = cx.tcx;
    let def = part.def;
    let msg = finding_message(tcx, def, sets, part.size, part.total);
    let site = site.map(|(site, args)| {
        (
            site.source_callsite(),
            format!(
                "one of the {sets} argument sets: {}",
                render_args(tcx, def, args)
            ),
        )
    });
    let rest = part_note(part, min_statements);
    let inline = inline_note(tcx, def);
    let help = format!(
        "move these statements into a separate non-generic function that takes what they read \
         and returns what they produce, and call it here. `{}` keeps its signature. They are \
         compiled once only if the compiler keeps the new function out of line",
        tcx.def_path_str(def)
    );
    let signature_help = part
        .signature_help
        .as_ref()
        .map(|line| format!("`{}` {line}", tcx.def_path_str(def)));
    emit_hir_then(
        cx,
        GENERIC_BODY_NOT_GENERIC,
        tcx.local_def_id_to_hir_id(def),
        tcx.def_span(def),
        msg,
        |diag| {
            match part.site {
                Some(SourceSpan::Whole(span)) => {
                    diag.span_note(span, format!("these statements {rest}"));
                }
                Some(SourceSpan::Start(span)) => {
                    let size = part.size;
                    diag.span_note(span, format!("the {size} statements starting here {rest}"));
                }
                None => {
                    diag.note(format!("these statements {rest}"));
                }
            }
            if let Some((site, note)) = site {
                diag.span_note(site, note);
            }
            if let Some(inline) = inline {
                diag.note(inline);
            }
            diag.help(help);
            if let Some(signature_help) = signature_help {
                diag.help(signature_help);
            }
        },
    );
}
