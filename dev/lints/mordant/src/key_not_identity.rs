use crate::baseline::emit;
use rustc_hir::def_id::{DefId, LOCAL_CRATE};
use rustc_hir::{Expr, ExprKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, Ty};
use rustc_span::sym;

use crate::MordantConfig;

rustc_session::declare_lint! {
    /// Flags a map keyed on something that does not identify the thing it
    /// names, so two different things can share a key and overwrite each
    /// other, or one thing can land under several: a type or method the
    /// project lists as not an identity, a float's `to_bits()`, a pointer
    /// cast to an integer. Which types and expression forms count is declared
    /// per project in `dylint.toml`. With no configuration the lint is
    /// silent.
    pub KEY_NOT_IDENTITY,
    Warn,
    "map keyed on a value that does not identify what it names"
}

/// A key-expression form `key-not-identity-forms` can opt into. These variants
/// are the whole of the accepted vocabulary, so a spelling this lint honours
/// exists in exactly one place; a name that is not one of them selects no form,
/// as it did when this was a pair of string comparisons.
#[derive(Clone, Copy, PartialEq, Eq)]
enum KeyForm {
    /// `f.to_bits()` where `f` is a float.
    ToBits,
    /// A raw pointer cast to an integer.
    PtrCast,
}

impl KeyForm {
    fn parse(name: &str) -> Option<Self> {
        match name {
            "to-bits" => Some(KeyForm::ToBits),
            "ptr-cast" => Some(KeyForm::PtrCast),
            _ => None,
        }
    }
}

pub struct KeyNotIdentity {
    config: &'static MordantConfig,
    /// The opted-in forms. Membership, not one flag per form: a form the
    /// project did not name is absent rather than false.
    forms: Vec<KeyForm>,
}

rustc_session::impl_lint_pass!(KeyNotIdentity => [KEY_NOT_IDENTITY]);

/// Methods whose receiver is a keyed collection and whose first argument (or
/// type parameter) is in key position.
const KEY_METHODS: &[&str] = &[
    "insert",
    "entry",
    "get",
    "contains_key",
    "contains",
    "remove",
];

impl KeyNotIdentity {
    pub fn new(config: &'static MordantConfig) -> Self {
        Self {
            config,
            forms: config
                .key_not_identity_forms
                .iter()
                .filter_map(|f| KeyForm::parse(f))
                .collect(),
        }
    }

    fn is_denied<'tcx>(&self, cx: &LateContext<'tcx>, ty: Ty<'tcx>) -> Option<String> {
        configured(
            cx,
            ty.peel_refs().ty_adt_def()?.did(),
            &self.config.key_not_identity_types,
        )
    }

    fn is_fixing<'tcx>(&self, cx: &LateContext<'tcx>, ty: Ty<'tcx>) -> bool {
        ty.peel_refs().ty_adt_def().is_some_and(|adt| {
            configured(cx, adt.did(), &self.config.key_not_identity_fixes).is_some()
        })
    }

    /// A direct hit, or (opt-in) a denied type inside a tuple or one level of
    /// struct fields with no declared identity-fixing component beside it:
    /// how the message introduces the key, and the denied type in it.
    fn denied_type_name<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        key_ty: Ty<'tcx>,
    ) -> Option<(String, String)> {
        if let Some(path) = self.is_denied(cx, key_ty) {
            return Some((
                format!("`{path}`, which this project's config says is not an identity"),
                path,
            ));
        }
        if !self.config.key_not_identity_composite {
            return None;
        }
        let components: Vec<Ty<'_>> = match key_ty.peel_refs().kind() {
            ty::Tuple(elems) => elems.iter().collect(),
            ty::Adt(adt, args) if adt.is_struct() => adt
                .non_enum_variant()
                .fields
                .iter()
                .map(|f| f.ty(cx.tcx, args).skip_normalization())
                .collect(),
            _ => return None,
        };
        let hit = components.iter().find_map(|t| self.is_denied(cx, *t))?;
        if components.iter().any(|t| self.is_fixing(cx, *t)) {
            return None;
        }
        Some((
            format!(
                "`{key_ty}`, which contains a `{hit}`, and this project's config says `{hit}` is not an identity"
            ),
            hit,
        ))
    }
}

/// `did`'s def path when `list` names it, either as that path or with the
/// local crate's name in front.
fn configured(cx: &LateContext<'_>, did: DefId, list: &[String]) -> Option<String> {
    let path = cx.tcx.def_path_str(did);
    let with_crate = did
        .is_local()
        .then(|| format!("{}::{path}", cx.tcx.crate_name(LOCAL_CRATE)));
    list.iter()
        .any(|d| *d == path || Some(d) == with_crate.as_ref())
        .then_some(path)
}

/// The key type of a keyed std/indexmap collection, or None.
fn keyed_collection_key<'tcx>(cx: &LateContext<'tcx>, recv_ty: Ty<'tcx>) -> Option<Ty<'tcx>> {
    let ty::Adt(adt, args) = recv_ty.peel_refs().kind() else {
        return None;
    };
    let did = adt.did();
    let is_std_keyed = [
        sym::HashMap,
        sym::HashSet,
        sym::BTreeMap,
        clippy_utils::sym::BTreeSet,
    ]
    .iter()
    .any(|s| cx.tcx.is_diagnostic_item(*s, did));
    let is_indexmap = cx.tcx.crate_name(did.krate).as_str() == "indexmap";
    if is_std_keyed || is_indexmap {
        args.types().next()
    } else {
        None
    }
}

fn is_float_to_bits(cx: &LateContext<'_>, expr: &Expr<'_>) -> bool {
    if let ExprKind::MethodCall(seg, recv, [], _) = expr.kind
        && seg.ident.as_str() == "to_bits"
    {
        cx.typeck_results()
            .expr_ty(recv)
            .peel_refs()
            .is_floating_point()
    } else {
        false
    }
}

fn is_ptr_to_int_cast(cx: &LateContext<'_>, expr: &Expr<'_>) -> bool {
    if let ExprKind::Cast(inner, _) = expr.kind {
        cx.typeck_results().expr_ty(inner).is_raw_ptr()
            && cx.typeck_results().expr_ty(expr).is_integral()
    } else {
        false
    }
}

impl<'tcx> LateLintPass<'tcx> for KeyNotIdentity {
    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let ExprKind::MethodCall(seg, recv, args, _) = expr.kind else {
            return;
        };
        if !KEY_METHODS.contains(&seg.ident.as_str()) {
            return;
        }
        let recv_ty = cx.typeck_results().expr_ty_adjusted(recv);
        let Some(key_ty) = keyed_collection_key(cx, recv_ty) else {
            return;
        };

        if let Some((shown, denied)) = self.denied_type_name(cx, key_ty) {
            let fixes: Vec<String> = self
                .config
                .key_not_identity_fixes
                .iter()
                .map(|f| format!("`{f}`"))
                .collect();
            let help = if fixes.is_empty() {
                format!("key on something that does identify the value instead of its `{denied}`")
            } else {
                format!(
                    "key on something that does identify the value. The config names {} as what makes a `{denied}` unambiguous, so add it to the key",
                    fixes.join(", "),
                )
            };
            emit(
                cx,
                KEY_NOT_IDENTITY,
                expr.span,
                format!(
                    "this map is keyed on {shown}. Two different things can share a `{denied}` and overwrite each other"
                ),
                help,
            );
            return;
        }

        // Expression-form checks apply where the key value itself is written:
        // `insert` and `entry`.
        if !matches!(seg.ident.as_str(), "insert" | "entry") {
            return;
        }
        let Some(key_expr) = args.first() else {
            return;
        };
        if self.forms.contains(&KeyForm::ToBits) && is_float_to_bits(cx, key_expr) {
            emit(
                cx,
                KEY_NOT_IDENTITY,
                key_expr.span,
                "this map is keyed on a float's `to_bits()`. For a NaN-boxed or interned value those bits are a pointer, so two equal values can land under different keys",
                "key on the value itself or its handle, not on its bit pattern",
            );
        }
        // Project-declared methods whose result is not an identity (e.g. a
        // NaN-boxing `Value::to_bits`, where boxed values yield pointer bits).
        let deny_methods = &self.config.key_not_identity_methods;
        if !deny_methods.is_empty()
            && let ExprKind::MethodCall(kseg, _, _, _) = key_expr.kind
            && let Some(mdid) = cx.typeck_results().type_dependent_def_id(key_expr.hir_id)
            && configured(cx, mdid, deny_methods).is_some()
        {
            let method = kseg.ident;
            emit(
                cx,
                KEY_NOT_IDENTITY,
                key_expr.span,
                format!(
                    "this map is keyed on the result of `{method}()`, which this project's config says is not an identity. Two equal things can land under different keys"
                ),
                format!("key on what identifies the value, not on what `{method}()` returns"),
            );
        }
        if self.forms.contains(&KeyForm::PtrCast) && is_ptr_to_int_cast(cx, key_expr) {
            emit(
                cx,
                KEY_NOT_IDENTITY,
                key_expr.span,
                "this map is keyed on a pointer cast to an integer, which names an allocation and not a value. Equal values at two addresses get two entries, and a freed and reused address looks like an old one",
                "key on the value, or on an id assigned when it is created, not on its address",
            );
        }
    }
}
