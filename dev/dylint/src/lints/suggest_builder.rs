use clippy_utils::diagnostics::span_lint_and_help;
use clippy_utils::ty::implements_trait;
use rustc_hir::{GenericParamKind, Item, ItemKind, LangItem, VariantData};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Symbol, sym};

use crate::config::SuggestBuilderConfig;

rustc_session::declare_lint! {
    /// Suggests using a `#[builder]` constructor in a `#[bon] impl` for structs with many named fields.
    pub SUGGEST_BUILDER,
    Warn,
    "suggests using a function builder for structs with many fields"
}

pub struct SuggestBuilder {
    threshold: usize,
    skip_derives: Vec<Symbol>,
}

impl SuggestBuilder {
    pub fn new() -> Self {
        let config: SuggestBuilderConfig = dylint_linting::config_or_default("suggest_builder");
        Self {
            threshold: config.threshold,
            skip_derives: config
                .skip_derives
                .iter()
                .map(|s| Symbol::intern(s))
                .collect(),
        }
    }
}

rustc_session::impl_lint_pass!(SuggestBuilder => [SUGGEST_BUILDER]);

impl<'tcx> LateLintPass<'tcx> for SuggestBuilder {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }
        let ItemKind::Struct(ident, generics, variant_data) = &item.kind else {
            return;
        };
        // Skip structs with lifetime parameters -- they represent borrowed
        // views, visitors, or traversal contexts where a builder is
        // structurally inappropriate.
        if generics
            .params
            .iter()
            .any(|p| matches!(p.kind, GenericParamKind::Lifetime { .. }))
        {
            return;
        }
        // Skip `#[repr(C)]` structs -- layout is dictated by FFI, a builder
        // is structurally inappropriate.
        let adt_def = cx.tcx.adt_def(item.owner_id);
        if adt_def.repr().c() {
            return;
        }
        let VariantData::Struct { fields, .. } = variant_data else {
            return;
        };
        // Don't count `PhantomData` fields (including variance markers like
        // `PhantomData<*const T>`, `PhantomData<fn(T)>`, etc.) -- they aren't
        // real from a construction-ergonomics standpoint.
        let field_count = fields
            .iter()
            .filter(|f| {
                let ty = cx.tcx.type_of(f.def_id).instantiate_identity().skip_norm_wip();
                !ty.ty_adt_def()
                    .is_some_and(|adt| cx.tcx.is_lang_item(adt.did(), LangItem::PhantomData))
            })
            .count();
        if field_count < self.threshold {
            return;
        }
        if super::has_bon_builder(ident.name) {
            return;
        }
        // Skip structs named `*Builder` -- they ARE builders, not builder
        // candidates.
        if ident.name.as_str().ends_with("Builder") {
            return;
        }
        // Skip structs that derive any trait in the configured `skip_derives`
        // list (default: Default, Queryable, Insertable, Selectable).
        if super::has_any_derive(ident.name, &self.skip_derives) {
            return;
        }
        // Skip structs that implement `Default` (derived or manual).
        let ty = cx.tcx.type_of(item.owner_id).instantiate_identity().skip_norm_wip();
        if let Some(default_id) = cx.tcx.get_diagnostic_item(rustc_span::sym::Default)
            && implements_trait(cx, ty, default_id, &[])
        {
            return;
        }
        // Skip structs with no inherent constructors (associated fns whose
        // return type contains Self). Without a constructor the struct is
        // built via literals -- a builder suggestion is noise.
        let struct_def_id = item.owner_id.to_def_id();
        let has_ctor = cx.tcx.inherent_impls(struct_def_id).iter().any(|impl_id| {
            cx.tcx
                .associated_items(*impl_id)
                .in_definition_order()
                .filter(|assoc| matches!(assoc.kind, ty::AssocKind::Fn { .. }))
                .any(|assoc| {
                    let ret_ty = cx
                        .tcx
                        .fn_sig(assoc.def_id)
                        .instantiate_identity()
                        .skip_norm_wip()
                        .output()
                        .skip_binder();
                    // Peel only the allowed constructor wrappers -- `Result<_, _>`
                    // (taking the `Ok` arg) and `Box<_>` (taking the inner arg) --
                    // then require an exact match to the struct's `Self` type.
                    // A return type that merely *contains* `Self` as a nested or
                    // generic arg (e.g. `Option<&Self>`) is NOT a constructor.
                    #[expect(
                        clippy::wildcard_enum_match_arm,
                        reason = "TyKind has many variants; only Result and Box are unwrapped"
                    )]
                    let inner_ty = match ret_ty.kind() {
                        // `Result<_, _>` (the `Ok` arg) and `Box<_>` (the inner
                        // arg) both carry the constructed type at generic arg 0.
                        ty::Adt(adt, args)
                            if cx.tcx.is_diagnostic_item(sym::Result, adt.did())
                                || cx.tcx.is_lang_item(adt.did(), LangItem::OwnedBox) =>
                        {
                            args.type_at(0)
                        }
                        _ => ret_ty,
                    };
                    inner_ty
                        .ty_adt_def()
                        .is_some_and(|adt| adt.did() == struct_def_id)
                })
        });
        if !has_ctor {
            return;
        }
        span_lint_and_help(
            cx,
            SUGGEST_BUILDER,
            item.span,
            format!(
                "struct `{ident}` has {field_count} fields; consider exposing a `#[builder]` constructor"
            ),
            None,
            "use the function builder style: `#[bon] impl` plus `#[builder] fn new(...) -> Self`",
        );
    }
}
