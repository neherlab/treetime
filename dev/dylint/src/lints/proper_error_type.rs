//! Flags unstructured error types in public API signatures.
//!
//! The TreeTime project standardizes on `eyre` / `color_eyre` for error
//! propagation: fallible functions return `Result<_, eyre::Report>` and build
//! errors with the `make_error!` / `make_report!` helper macros. This lint flags
//! public and `pub(crate)` functions that instead return a stringly-typed error
//! (`String`, `&str`, `Cow<str>`), an untyped `Box<dyn Error>`, or a competing
//! erased error crate (`anyhow::Error`, `miette::Report`). `eyre::Report` itself
//! is the approved type and is never flagged.

use clippy_utils::diagnostics::span_lint_and_help;
use clippy_utils::{is_def_id_trait_method, is_entrypoint_fn, is_in_cfg_test, return_ty};
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, FnDecl, LangItem};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, ExistentialPredicate, Ty};
use rustc_span::def_id::CRATE_DEF_ID;
use rustc_span::{Symbol, sym};

const HELP: &str =
    "return `eyre::Report` (or `color_eyre::Report`) and construct errors with `make_error!` / `make_report!`";

enum UnstructuredKind {
    Basic(&'static str),
    ErasedCrate {
        crate_name: &'static str,
        type_name: &'static str,
    },
}

rustc_session::declare_lint! {
    /// Flags error types in public APIs that are stringly-typed, untyped, or
    /// from a competing erased error crate instead of the project's `eyre`.
    pub PROPER_ERROR_TYPE,
    Warn,
    "flags improper error types in public API signatures"
}

pub struct ProperErrorType {
    sym_error: Symbol,
    sym_anyhow: Symbol,
    sym_miette: Symbol,
    sym_miette_report: Symbol,
}

impl Default for ProperErrorType {
    fn default() -> Self {
        Self {
            sym_error: Symbol::intern("Error"),
            sym_anyhow: Symbol::intern("anyhow"),
            sym_miette: Symbol::intern("miette"),
            sym_miette_report: Symbol::intern("Report"),
        }
    }
}

impl ProperErrorType {
    /// Returns the visibility when an item has explicit `pub` or `pub(crate)` visibility.
    ///
    /// Private items at the crate root share the same semantic visibility
    /// (`Restricted(CRATE_DEF_ID)`) as `pub(crate)` items, so we require the HIR
    /// node to carry a non-empty `vis_span` (i.e. the user actually wrote a
    /// visibility keyword).
    fn vis_if_at_least_pub_crate(
        cx: &LateContext<'_>,
        def_id: rustc_hir::def_id::DefId,
    ) -> Option<ty::Visibility<rustc_hir::def_id::DefId>> {
        let vis = cx.tcx.visibility(def_id);
        if vis.is_public() {
            return Some(vis);
        }
        if vis != ty::Visibility::Restricted(CRATE_DEF_ID.to_def_id()) {
            return None;
        }
        // Semantic visibility is pub(crate); verify an explicit keyword exists.
        #[expect(
            clippy::wildcard_enum_match_arm,
            reason = "only Item and ImplItem carry a visibility span; all other HIR nodes \
                      are definitions whose visibility is inherited from their parent"
        )]
        def_id
            .as_local()
            .is_some_and(|local| match cx.tcx.hir_node_by_def_id(local) {
                rustc_hir::Node::Item(item) => !item.vis_span.is_empty(),
                rustc_hir::Node::ImplItem(ii) => ii.vis_span().is_some(),
                _ => false,
            })
            .then_some(vis)
    }

    fn classify_unstructured<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        err_ty: Ty<'tcx>,
    ) -> Option<UnstructuredKind> {
        #[expect(
            clippy::wildcard_enum_match_arm,
            reason = "ty::TyKind has too many variants to enumerate; all unrecognised kinds correctly return None"
        )]
        match err_ty.kind() {
            ty::Ref(_, inner, _) if inner.is_str() => Some(UnstructuredKind::Basic("&str")),
            ty::Adt(adt, args) => {
                let did = adt.did();
                if cx.tcx.is_lang_item(did, LangItem::String) {
                    return Some(UnstructuredKind::Basic("String"));
                }
                if cx.tcx.is_diagnostic_item(sym::Cow, did) && args.type_at(1).is_str() {
                    return Some(UnstructuredKind::Basic("Cow<'_, str>"));
                }
                if cx.tcx.is_lang_item(did, LangItem::OwnedBox)
                    && let ty::Dynamic(preds, ..) = args.type_at(0).kind()
                    && let Some(error_trait_id) = cx.tcx.get_diagnostic_item(self.sym_error)
                    && preds.iter().any(|pred| {
                        matches!(pred.skip_binder(), ExistentialPredicate::Trait(t) if t.def_id == error_trait_id)
                    })
                {
                    return Some(UnstructuredKind::Basic("Box<dyn Error>"));
                }
                let crate_name = cx.tcx.crate_name(did.krate);
                let item_name = cx.tcx.item_name(did);
                if crate_name == self.sym_anyhow && item_name == self.sym_error {
                    Some(UnstructuredKind::ErasedCrate { crate_name: "anyhow", type_name: "Error" })
                } else if crate_name == self.sym_miette && item_name == self.sym_miette_report {
                    Some(UnstructuredKind::ErasedCrate { crate_name: "miette", type_name: "Report" })
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

rustc_session::impl_lint_pass!(ProperErrorType => [PROPER_ERROR_TYPE]);

impl<'tcx> LateLintPass<'tcx> for ProperErrorType {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        decl: &'tcx FnDecl<'tcx>,
        _body: &'tcx Body<'tcx>,
        span: rustc_span::Span,
        def_id: rustc_hir::def_id::LocalDefId,
    ) {
        if span.from_expansion()
            || is_def_id_trait_method(cx, def_id)
            || is_entrypoint_fn(cx, def_id.to_def_id())
            || is_in_cfg_test(cx.tcx, cx.tcx.local_def_id_to_hir_id(def_id))
            || matches!(kind, FnKind::Closure)
        {
            return;
        }

        let Some(vis) = Self::vis_if_at_least_pub_crate(cx, def_id.to_def_id()) else {
            return;
        };

        let ret_ty = return_ty(cx, rustc_hir::OwnerId { def_id });
        let ty::Adt(adt, args) = ret_ty.kind() else {
            return;
        };
        if !cx.tcx.is_diagnostic_item(sym::Result, adt.did()) {
            return;
        }
        let err_ty = args.type_at(1);

        let Some(kind) = self.classify_unstructured(cx, err_ty) else {
            return;
        };

        let ret_span = decl.output.span();

        match kind {
            UnstructuredKind::Basic(name) => {
                let vis_label = if vis.is_public() {
                    "public"
                } else {
                    "`pub(crate)`"
                };
                span_lint_and_help(
                    cx,
                    PROPER_ERROR_TYPE,
                    ret_span,
                    format!(
                        "{vis_label} function returns `Result<_, {name}>` -- use `eyre::Report`"
                    ),
                    None,
                    HELP,
                );
            }
            UnstructuredKind::ErasedCrate { crate_name, type_name } => {
                if !cx.tcx.effective_visibilities(()).is_reachable(def_id) {
                    return;
                }
                span_lint_and_help(
                    cx,
                    PROPER_ERROR_TYPE,
                    ret_span,
                    format!(
                        "effectively public function returns `{crate_name}::{type_name}` -- the project standardizes on `eyre`"
                    ),
                    None,
                    HELP,
                );
            }
        }
    }
}
