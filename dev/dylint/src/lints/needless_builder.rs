use clippy_utils::diagnostics::span_lint_and_help;
use rustc_hir::{Item, ItemKind, LangItem, VariantData};
use rustc_lint::{LateContext, LateLintPass};

use crate::config::NeedlessBuilderConfig;

rustc_session::declare_lint! {
    /// Warns when `bon::Builder` is derived on a struct with very few fields.
    pub NEEDLESS_BUILDER,
    Warn,
    "warns when `bon::Builder` is used on structs with few fields"
}

pub struct NeedlessBuilder {
    threshold: usize,
}

impl NeedlessBuilder {
    pub fn new() -> Self {
        let config: NeedlessBuilderConfig = dylint_linting::config_or_default("needless_builder");
        Self {
            threshold: config.threshold,
        }
    }
}

rustc_session::impl_lint_pass!(NeedlessBuilder => [NEEDLESS_BUILDER]);

impl<'tcx> LateLintPass<'tcx> for NeedlessBuilder {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }
        let ItemKind::Struct(ident, _, variant_data) = &item.kind else {
            return;
        };
        let VariantData::Struct { fields, .. } = variant_data else {
            return;
        };
        // Cheap thread-local gate first: only structs deriving `bon::Builder`
        // can be flagged, so don't pay for the per-field `type_of` queries below
        // on every other struct in the crate.
        if !super::has_bon_builder(ident.name) {
            return;
        }
        // Don't count `PhantomData` markers as real fields -- a struct of two
        // real fields plus a `PhantomData<T>` should still be flagged.
        let field_count = fields
            .iter()
            .filter(|f| {
                let ty = cx.tcx.type_of(f.def_id).instantiate_identity().skip_norm_wip();
                !ty.ty_adt_def()
                    .is_some_and(|adt| cx.tcx.is_lang_item(adt.did(), LangItem::PhantomData))
            })
            .count();
        if field_count > self.threshold {
            return;
        }
        span_lint_and_help(
            cx,
            NEEDLESS_BUILDER,
            item.span,
            format!(
                "struct `{ident}` has only {field_count} fields; `bon::Builder` may be unnecessary",
            ),
            None,
            "consider using a plain constructor or struct literal instead",
        );
    }
}
