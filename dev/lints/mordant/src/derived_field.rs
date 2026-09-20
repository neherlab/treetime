use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, StructTailExpr};
use rustc_lint::{LateContext, LateLintPass};
use rustc_span::{Span, Symbol};

use crate::MordantConfig;
use crate::adt_facts::{has_fixed_repr, has_positional_fields};
use crate::baseline::emit_with_note;
use crate::enum_facts::ctor_literal_variant;
use crate::hir_shapes::{assigned_adt_field, peel_blocks_unsafe};

rustc_session::declare_lint! {
    /// Flags a field that always has the same value for a given value of a
    /// sibling field, in every place the type is built. It is a stored copy
    /// of something the sibling decides, and nothing rejects a mismatched
    /// pair written by hand. A method on the sibling's enum replaces it.
    ///
    /// Fires only when one of the two is a variant of an enum *this crate
    /// defined* at every site — a closed set of states, which is what makes
    /// the repair a method on that enum. Literal and named-constant columns
    /// biject whenever a small table's rows differ, which is a property of
    /// the table and not of the type; see `Val::decides`.
    ///
    /// A field that is also assigned somewhere (`x.f = …`) is skipped: the
    /// constructors are then not the only thing deciding it, and a value the
    /// literals always pair one way may be re-paired later.
    ///
    /// Silent on any type with an explicit `repr`, on foreign types, and
    /// below `derived-field-min-sites` construction sites.
    pub DERIVED_FIELD,
    Warn,
    "a field whose value is decided by a sibling field"
}

pub struct DerivedField {
    min_sites: usize,
    seen: HashMap<DefId, Vec<Site>>,
    /// (variant, field) pairs written by an assignment rather than a literal.
    assigned: HashMap<DefId, BTreeSet<Symbol>>,
}

rustc_session::impl_lint_pass!(DerivedField => [DERIVED_FIELD]);

/// One site cannot exhibit a correspondence and two is the least that can.
const FLOOR: usize = 2;

impl DerivedField {
    pub fn new(config: &MordantConfig) -> Self {
        Self {
            min_sites: config.derived_field_min_sites.max(FLOOR),
            seen: HashMap::new(),
            assigned: HashMap::new(),
        }
    }
}

struct Site {
    span: Span,
    fields: BTreeMap<Symbol, Val>,
}

/// What a field initialiser evaluates to, as far as a name can be put on it.
///
/// A `Variant` records *which* variant and never its payload: `Some(n)` and
/// `Some(m)` are the same fact here, which is what makes an `Option` beside an
/// enum legible as one correspondence rather than as many.
#[derive(Clone, PartialEq, Eq, Hash)]
enum Val {
    Variant(DefId),
    NamedConst(DefId),
    Lit(String),
}

impl Val {
    /// Can this value *decide* a sibling — is it drawn from a closed set this
    /// crate defined?
    ///
    /// Only a local enum variant is. That is a stronger test than "has a
    /// name", and each of the three things it rules out was measured firing
    /// before it was added:
    ///
    /// * **Bare literals.** Two literal columns are in bijection whenever
    ///   their rows happen to differ, which is a property of a table with few
    ///   rows and not of the type. Every such pairing on the tree this lint
    ///   was measured against was a fixture.
    /// * **Named constants.** A two-row table pairing distinct strings with a
    ///   per-row constant satisfies "one side is named" for free, and the
    ///   advice to replace it with a method is wrong there — the constant is
    ///   the row's own datum, not a restatement of its neighbour.
    /// * **Foreign variants**, which is what keeps this off
    ///   `options_as_enum`' ground: two `Option` fields set to opposite
    ///   presences across two constructors biject, but `Some` and `None` are
    ///   `core`'s names and that shape already has a lint.
    ///
    /// What is left is the case the fix is actually written for: an enum
    /// names a closed set of states, a sibling restates one of them, and the
    /// repair is a method on the enum — `Ceiling::limit()`, or folding the
    /// field into the variant as `Wide { remaining }`.
    fn decides(&self) -> bool {
        matches!(self, Val::Variant(d) if d.is_local())
    }
}

/// Strip the wrappers that do not change which value an expression names.
fn peel<'tcx>(e: &'tcx Expr<'tcx>) -> &'tcx Expr<'tcx> {
    let e = peel_blocks_unsafe(e);
    match e.kind {
        ExprKind::AddrOf(_, _, inner) => peel(inner),
        _ => e,
    }
}

fn classify<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Option<Val> {
    let e = peel(e);
    // `Some(x)`, `Wide(n)` — the variant is the fact, the payload is not.
    if let Some(v) = ctor_literal_variant(cx, e) {
        return Some(Val::Variant(v));
    }
    match e.kind {
        ExprKind::Path(ref qpath) => match cx.qpath_res(qpath, e.hir_id) {
            Res::Def(DefKind::Const { .. } | DefKind::AssocConst { .. }, did) => {
                Some(Val::NamedConst(did))
            }
            _ => None,
        },
        ExprKind::Lit(lit) => Some(Val::Lit(format!("{:?}", lit.node))),
        _ => None,
    }
}

impl DerivedField {
    /// `s.f = …`, through any number of derefs and autoderefs: `f` of `s`'s
    /// struct is decided by more than its literals.
    fn note_assignment<'tcx>(&mut self, cx: &LateContext<'tcx>, place: &'tcx Expr<'tcx>) {
        let Some((adt, field, _)) = assigned_adt_field(cx, place) else {
            return;
        };
        if !adt.is_struct() || !adt.did().is_local() {
            return;
        }
        self.assigned
            .entry(adt.non_enum_variant().def_id)
            .or_default()
            .insert(field.name);
    }
}

impl<'tcx> LateLintPass<'tcx> for DerivedField {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        if let ExprKind::Assign(place, ..) | ExprKind::AssignOp(_, place, _) = expr.kind {
            self.note_assignment(cx, place);
            return;
        }
        // A `..base` literal overrides part of a value that already exists, so
        // the fields written are the ones deliberately being made to differ.
        let ExprKind::Struct(qpath, fields, StructTailExpr::None) = expr.kind else {
            return;
        };
        // Typeck rather than the path, so an alias or `Self` resolves.
        let Some(adt) = cx.typeck_results().expr_ty(expr).ty_adt_def() else {
            return;
        };
        if !adt.did().is_local() {
            return;
        }
        let variant = adt.variant_of_res(cx.qpath_res(qpath, expr.hir_id));
        // A wire record legitimately restates what a sibling implies.
        if has_fixed_repr(adt) || has_positional_fields(variant) {
            return;
        }
        let mut vals = BTreeMap::new();
        for f in fields {
            if let Some(v) = classify(cx, f.expr) {
                vals.insert(f.ident.name, v);
            }
        }
        self.seen.entry(variant.def_id).or_default().push(Site {
            span: expr.span,
            fields: vals,
        });
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let mut keys: Vec<DefId> = self
            .seen
            .iter()
            .filter(|(_, s)| s.len() >= self.min_sites)
            .map(|(d, _)| *d)
            .collect();
        // Deterministic order: findings are compared against a baseline.
        keys.sort_by_key(|d| cx.tcx.def_span(*d));

        for did in keys {
            let sites = &self.seen[&did];
            // Fields every site gave a name to. A site that left one unknown
            // cannot witness a correspondence, so the whole field is dropped
            // rather than the pairing being read off a subset.
            let mut common: BTreeSet<Symbol> = sites[0].fields.keys().copied().collect();
            for s in &sites[1..] {
                common.retain(|k| s.fields.contains_key(k));
            }
            if let Some(assigned) = self.assigned.get(&did) {
                common.retain(|f| !assigned.contains(f));
            }
            let common: Vec<Symbol> = common.into_iter().collect();
            for i in 0..common.len() {
                for j in (i + 1)..common.len() {
                    let (a, b) = (common[i], common[j]);
                    let va: Vec<&Val> = sites.iter().map(|s| &s.fields[&a]).collect();
                    let vb: Vec<&Val> = sites.iter().map(|s| &s.fields[&b]).collect();
                    if !(va.iter().all(|v| v.decides()) || vb.iter().all(|v| v.decides())) {
                        continue;
                    }
                    let da: HashSet<&&Val> = va.iter().collect();
                    let db: HashSet<&&Val> = vb.iter().collect();
                    // One distinct value on either side is a constant, not a
                    // correspondence.
                    if da.len() < 2 || db.len() < 2 {
                        continue;
                    }
                    if !bijective(&va, &vb) {
                        continue;
                    }
                    // The side drawn from a local enum decides the other; the
                    // repair is a method on that enum.
                    let (deciding, derived, by) = if va.iter().all(|v| v.decides()) {
                        (a, b, va[0])
                    } else {
                        (b, a, vb[0])
                    };
                    let Val::Variant(variant) = by else { continue };
                    let enum_name = cx.tcx.item_name(cx.tcx.parent(*variant));
                    emit_with_note(
                        cx,
                        DERIVED_FIELD,
                        cx.tcx.def_span(did),
                        format!(
                            "`{derived}` always has the same value for a given `{deciding}`, in all \
                             {} places `{}` is built. It is a stored copy of something `{deciding}` \
                             decides",
                            sites.len(),
                            cx.tcx.def_path_str(did),
                        ),
                        sites[0].span,
                        "one of the constructions",
                        format!(
                            "add `fn {derived}(&self)` to `{enum_name}` and delete the `{derived}` field"
                        ),
                    );
                }
            }
        }
    }
}

/// Do `a` and `b` agree one-for-one — each value of one always beside the same
/// value of the other, in both directions?
fn bijective(a: &[&Val], b: &[&Val]) -> bool {
    let mut fwd: HashMap<&Val, &Val> = HashMap::new();
    let mut rev: HashMap<&Val, &Val> = HashMap::new();
    for (x, y) in a.iter().zip(b.iter()) {
        if *fwd.entry(x).or_insert(y) != *y {
            return false;
        }
        if *rev.entry(y).or_insert(x) != *x {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use rustc_hir::def_id::{CrateNum, DefId, DefIndex};

    use super::{Val, bijective};

    fn lit(s: &str) -> Val {
        Val::Lit(s.to_string())
    }

    fn def(krate: u32, index: u32) -> DefId {
        DefId {
            krate: CrateNum::from_u32(krate),
            index: DefIndex::from_u32(index),
        }
    }

    #[test]
    fn one_for_one_is_bijective() {
        let (a, b) = (vec![lit("x"), lit("y")], vec![lit("1"), lit("2")]);
        let (ar, br): (Vec<&Val>, Vec<&Val>) = (a.iter().collect(), b.iter().collect());
        assert!(bijective(&ar, &br));
    }

    /// The same left value beside two different right values is the whole
    /// point: the pairing is not a function and the field is not a projection.
    #[test]
    fn a_repeated_left_with_two_rights_is_not() {
        let (a, b) = (vec![lit("x"), lit("x")], vec![lit("1"), lit("2")]);
        let (ar, br): (Vec<&Val>, Vec<&Val>) = (a.iter().collect(), b.iter().collect());
        assert!(!bijective(&ar, &br));
    }

    /// And the reverse direction, which a one-way functional check would pass.
    #[test]
    fn two_lefts_sharing_one_right_is_not() {
        let (a, b) = (vec![lit("x"), lit("y")], vec![lit("1"), lit("1")]);
        let (ar, br): (Vec<&Val>, Vec<&Val>) = (a.iter().collect(), b.iter().collect());
        assert!(!bijective(&ar, &br));
    }

    /// A bare literal decides nothing: a column of them is in bijection with
    /// its neighbour whenever a short table's rows differ.
    #[test]
    fn bare_literals_decide_nothing() {
        assert!(!lit("1").decides());
    }

    /// Nor does a constant this crate named, which is the two-row-table false
    /// positive this lint was measured emitting before `decides` narrowed.
    #[test]
    fn a_named_constant_decides_nothing() {
        assert!(!Val::NamedConst(def(0, 7)).decides());
    }

    /// A foreign variant is `Some`/`None`, which is `options_as_enum`' shape
    /// rather than this one.
    #[test]
    fn a_foreign_variant_decides_nothing() {
        assert!(!Val::Variant(def(2, 7)).decides());
    }

    #[test]
    fn a_local_enum_variant_decides() {
        assert!(Val::Variant(def(0, 7)).decides());
    }
}
