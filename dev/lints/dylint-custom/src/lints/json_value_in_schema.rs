//! Flags untyped JSON (`serde_json::Value`, `serde_json::Map`, or a container
//! of them) in the fields of a type that implements `schemars::JsonSchema`.
//! Such a field becomes an open schema, which the generated TypeScript types
//! render as `unknown`, so every client has to re-parse it by hand.

use clippy_utils::diagnostics::span_lint_and_help;
use rustc_hir::def::Res;
use rustc_hir::{Attribute, Item, ItemKind};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_middle::ty::{self, Ty};
use rustc_span::Span;
use rustc_span::def_id::DefId;

rustc_session::declare_lint! {
    /// Flags a field of a `JsonSchema` type whose type is `serde_json::Value`,
    /// `serde_json::Map`, or a container of them, and a `#[schemars(with = ...)]`
    /// that names one of these types. Use a typed struct or enum, or the named
    /// JSON type `JsonValue` when the data is free-form by design.
    pub JSON_VALUE_IN_SCHEMA,
    Warn,
    "untyped JSON in a field of a `JsonSchema` type -- the generated client sees `unknown`"
}

const JSON_SCHEMA_TRAIT: &str = "schemars::JsonSchema";

const UNTYPED_JSON: &[&str] = &["serde_json::Value", "serde_json::Map"];

const ALLOWED_TYPES: &[&str] = &["JsonValue", "SparseConfig", "AuspiceDocument"];

pub struct JsonValueInSchema;

impl JsonValueInSchema {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(JsonValueInSchema => [JSON_VALUE_IN_SCHEMA]);

impl<'tcx> LateLintPass<'tcx> for JsonValueInSchema {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        let ItemKind::Impl(impl_block) = &item.kind else {
            return;
        };
        let Some(of_trait) = &impl_block.of_trait else {
            return;
        };
        let Res::Def(_, trait_id) = of_trait.trait_ref.path.res else {
            return;
        };
        if !cx.tcx.def_path_str(trait_id).ends_with(JSON_SCHEMA_TRAIT) {
            return;
        }
        let self_ty = cx.tcx.type_of(item.owner_id).instantiate_identity().skip_norm_wip();
        let ty::Adt(adt, _) = self_ty.kind() else {
            return;
        };
        if !adt.did().is_local() {
            return;
        }
        let type_name = cx.tcx.item_name(adt.did());
        if ALLOWED_TYPES.contains(&type_name.as_str()) {
            return;
        }
        for field in adt.all_fields() {
            let field_ty = cx.tcx.type_of(field.did).instantiate_identity().skip_norm_wip();
            let span = cx.tcx.def_span(field.did);
            if is_untyped_json(cx, field_ty) || names_untyped_json_in_attribute(cx, field.did) {
                report(cx, span, field.name.as_str(), type_name.as_str());
            }
        }
    }
}

fn is_untyped_json<'tcx>(cx: &LateContext<'tcx>, ty: Ty<'tcx>) -> bool {
    ty.walk().any(|arg| {
        arg.as_type().is_some_and(|inner| {
            matches!(inner.kind(), ty::Adt(adt, _) if UNTYPED_JSON.contains(&cx.tcx.def_path_str(adt.did()).as_str()))
        })
    })
}

fn names_untyped_json_in_attribute(cx: &LateContext<'_>, field: DefId) -> bool {
    let Some(local) = field.as_local() else {
        return false;
    };
    let hir_id = cx.tcx.local_def_id_to_hir_id(local);
    cx.tcx.hir_attrs(hir_id).iter().any(|attr| {
        if !matches!(attr, Attribute::Unparsed(_)) {
            return false;
        }
        let Ok(snippet) = cx.sess().source_map().span_to_snippet(attr.span()) else {
            return false;
        };
        let compact: String = snippet.chars().filter(|c| !c.is_whitespace()).collect();
        compact.contains("schemars(")
            && compact.contains("with=\"")
            && (compact.contains("Value\"") || compact.contains("Map<"))
    })
}

fn report(cx: &LateContext<'_>, span: Span, field: &str, type_name: &str) {
    span_lint_and_help(
        cx,
        JSON_VALUE_IN_SCHEMA,
        span,
        format!("field `{field}` of `{type_name}` is untyped JSON in a schema type"),
        None,
        "give the field a typed struct or enum, or the named JSON type `JsonValue` when the data is free-form by design",
    );
}
