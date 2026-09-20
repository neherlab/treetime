use std::collections::{HashMap, HashSet};

use rustc_hir::def_id::DefId;
use rustc_hir::{BinOpKind, Expr, ExprKind, Pat};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;

use crate::adt_facts::inside_own_trait_impl;
use crate::baseline::emit_with_note;
use crate::enum_facts::{arm_variant, private_enum_of};

rustc_session::declare_lint! {
    /// Flags a crate-private enum variant that is constructed somewhere but
    /// no pattern anywhere singles it out (its own trait impls aside), so
    /// what it carries is never read and it is only ever caught by a
    /// catch-all. Trait
    /// impls (`Display`, `Debug`, `From`, derives) must match every variant to
    /// exist, so they prove nothing. A pattern anywhere else, including the
    /// enum's own inherent methods, is the crate reading the structure. An
    /// `==` / `!=` against a variant names it as a pattern would, and an `as`
    /// cast of the enum reads every variant, since the discriminant is all the
    /// structure a fieldless enum has. A fieldless variant that is the only
    /// unnamed one in its enum is reached by elimination and stays silent. A
    /// variant that fails all of this only ever reaches anyone through a
    /// catch-all shared with a sibling, or a string rendering.
    pub UNREAD_ERROR_VARIANT,
    Warn,
    "enum variant constructed but never distinguished by a pattern"
}

#[derive(Default)]
struct EnumFacts {
    /// variant -> first construction site.
    constructed: HashMap<DefId, rustc_span::Span>,
    /// Variants named by a pattern outside the enum's own impls.
    named: HashSet<DefId>,
}

#[derive(Default)]
pub struct UnreadErrorVariant {
    enums: HashMap<DefId, EnumFacts>,
}

rustc_session::impl_lint_pass!(UnreadErrorVariant => [UNREAD_ERROR_VARIANT]);

/// The private-enum variant `expr` constructs, if it is a variant path
/// (including a bare tuple constructor passed as a function), call or struct
/// expression, as `(enum, variant)`.
fn constructed_variant(cx: &LateContext<'_>, expr: &Expr<'_>) -> Option<(DefId, DefId)> {
    let variant = crate::enum_facts::constructed_variant(cx, expr)?;
    Some((private_enum_of(cx, variant)?, variant))
}

impl<'tcx> LateLintPass<'tcx> for UnreadErrorVariant {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        match expr.kind {
            // `x == E::V` singles out `V` exactly as `matches!(x, E::V)`
            // would. The derived `PartialEq` that runs is a trait impl, but
            // this comparison site is not, so it counts wherever a pattern
            // would.
            ExprKind::Binary(op, lhs, rhs) if matches!(op.node, BinOpKind::Eq | BinOpKind::Ne) => {
                for side in [lhs, rhs] {
                    if let Some((enum_did, variant)) =
                        constructed_variant(cx, clippy_utils::peel_ref_operators(cx, side))
                        && !inside_own_trait_impl(cx, expr.hir_id, enum_did)
                    {
                        self.enums
                            .entry(enum_did)
                            .or_default()
                            .named
                            .insert(variant);
                    }
                }
            }
            // `e as u8` observes the discriminant of whatever `e` holds, so
            // every variant of the enum is read by it.
            ExprKind::Cast(inner, _) => {
                if let ty::Adt(adt, _) = cx.typeck_results().expr_ty(inner).kind()
                    && adt.is_enum()
                    && let Some(local) = adt.did().as_local()
                    && !cx.effective_visibilities.is_exported(local)
                {
                    self.enums
                        .entry(adt.did())
                        .or_default()
                        .named
                        .extend(adt.variants().iter().map(|v| v.def_id));
                }
            }
            _ => {}
        }
        let Some((enum_did, variant)) = constructed_variant(cx, expr) else {
            return;
        };
        self.enums
            .entry(enum_did)
            .or_default()
            .constructed
            .entry(variant)
            .or_insert(expr.span);
    }

    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let Some(variant) = arm_variant(cx, pat) else {
            return;
        };
        let Some(enum_did) = private_enum_of(cx, variant) else {
            return;
        };
        if inside_own_trait_impl(cx, pat.hir_id, enum_did) {
            return;
        }
        self.enums
            .entry(enum_did)
            .or_default()
            .named
            .insert(variant);
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut unread = Vec::new();
        for (enum_did, facts) in &self.enums {
            // The crate must distinguish at least one variant by a counting
            // pattern; otherwise matching is simply not how this enum is
            // consumed, and "never matched" proves nothing about any variant.
            if facts.named.is_empty() {
                continue;
            }
            let variants = cx.tcx.adt_def(*enum_did).variants();
            let unnamed: Vec<_> = variants
                .iter()
                .filter(|v| !facts.named.contains(&v.def_id))
                .collect();
            // When every other variant is named, the last one is singled out
            // by elimination (`if s != Active && s != Inactive` reaches
            // exactly `Done`). That reads a fieldless variant completely; only
            // a payload could still go unread.
            if let [only] = unnamed.as_slice()
                && only.fields.is_empty()
            {
                continue;
            }
            for (variant, span) in &facts.constructed {
                if !facts.named.contains(variant) {
                    unread.push((*span, *variant, *enum_did));
                }
            }
        }
        unread.sort_by_key(|(span, ..)| span.lo());
        for (span, variant, enum_did) in unread {
            let path = cx.tcx.def_path_str(variant);
            let owner = cx.tcx.item_name(enum_did);
            emit_with_note(
                cx,
                UNREAD_ERROR_VARIANT,
                span,
                format!(
                    "`{path}` is constructed here, but no pattern anywhere singles it out (its own trait impls aside). What it carries is never read, and it is only ever caught by a catch-all"
                ),
                cx.tcx.def_span(variant),
                "the variant",
                format!(
                    "add an arm for `{path}` where `{owner}` is matched, or fold it into the variant it is lumped with"
                ),
            );
        }
    }
}
