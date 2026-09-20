use crate::adt_facts::result_err_ty;
use crate::baseline::emit;
use clippy_utils::ty::ty_from_hir_ty;
use rustc_hir::def_id::LocalDefId;
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, FnDecl, FnRetTy, TraitFn, TraitItem, TraitItemKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty::{self, Ty};
use rustc_span::{Span, Symbol, sym};

use crate::MordantConfig;

rustc_session::declare_lint! {
    /// Flags a public function whose error type is `String` (or `&str`,
    /// `Cow<str>`): callers can only tell failures apart by reading the
    /// message. An error enum gives them variants to match on.
    pub STRINGLY_ERROR,
    Warn,
    "public signature returns a Result with a string error type"
}

pub struct StringlyError {
    pub config: &'static MordantConfig,
}

rustc_session::impl_lint_pass!(StringlyError => [STRINGLY_ERROR]);

impl StringlyError {
    fn check_sig<'tcx>(&self, cx: &LateContext<'tcx>, def_id: LocalDefId, decl: &FnDecl<'tcx>) {
        if !cx.effective_visibilities.is_exported(def_id) {
            return;
        }
        let FnRetTy::Return(ret_hir_ty) = decl.output else {
            return;
        };
        let output = ty_from_hir_ty(cx, ret_hir_ty);
        let Some(err_ty) = result_err_ty(cx.tcx, output) else {
            return;
        };
        if let Some(desc) = self.stringy_desc(cx, err_ty) {
            let name = cx.tcx.item_name(def_id.to_def_id());
            emit(
                cx,
                STRINGLY_ERROR,
                ret_hir_ty.span,
                format!(
                    "`{name}` is public and returns `Result<_, {desc}>`. Callers can only tell failures apart by reading the message"
                ),
                "return an error enum with a variant per failure, and keep the text in its `Display` impl",
            );
        }
    }

    fn stringy_desc(&self, cx: &LateContext<'_>, err_ty: Ty<'_>) -> Option<&'static str> {
        match err_ty.kind() {
            ty::Adt(adt, args) => {
                if cx.tcx.is_lang_item(adt.did(), rustc_hir::LangItem::String) {
                    Some("String")
                } else if cx.tcx.is_diagnostic_item(sym::Cow, adt.did()) && args.type_at(1).is_str()
                {
                    Some("Cow<str>")
                } else if self.config.stringly_error_include_box_dyn
                    && adt.is_box()
                    && is_dyn_error(cx, args.type_at(0))
                {
                    Some("Box<dyn Error>")
                } else {
                    None
                }
            }
            ty::Ref(_, inner, _) if inner.is_str() => Some("&str"),
            _ => None,
        }
    }
}

fn is_dyn_error(cx: &LateContext<'_>, ty: Ty<'_>) -> bool {
    if let ty::Dynamic(preds, ..) = ty.kind()
        && let Some(principal) = preds.principal_def_id()
    {
        cx.tcx.get_diagnostic_item(Symbol::intern("Error")) == Some(principal)
    } else {
        false
    }
}

impl<'tcx> LateLintPass<'tcx> for StringlyError {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        decl: &'tcx FnDecl<'tcx>,
        _body: &'tcx Body<'tcx>,
        _span: Span,
        def_id: LocalDefId,
    ) {
        if matches!(kind, FnKind::Closure) {
            return;
        }
        // `fn main() -> Result<(), String>` in examples and binaries is not an API.
        if cx
            .tcx
            .entry_fn(())
            .is_some_and(|(did, _)| did == def_id.to_def_id())
        {
            return;
        }
        // A trait-impl method's signature is dictated by the trait. The finding,
        // if any, belongs on the trait definition.
        if let Some(assoc) = cx.tcx.opt_associated_item(def_id.to_def_id())
            && assoc.trait_item_def_id().is_some()
        {
            return;
        }
        self.check_sig(cx, def_id, decl);
    }

    fn check_trait_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx TraitItem<'tcx>) {
        // Required methods have no body, so check_fn never sees them.
        if let TraitItemKind::Fn(sig, TraitFn::Required(_)) = &item.kind {
            self.check_sig(cx, item.owner_id.def_id, sig.decl);
        }
    }
}
