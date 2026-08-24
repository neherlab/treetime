use clippy_utils::diagnostics::span_lint_and_help;
use clippy_utils::{is_def_id_trait_method, return_ty};
use rustc_hir::intravisit::FnKind;
use rustc_hir::{Body, FnDecl, Item, ItemKind};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::{Span, sym};

rustc_session::declare_lint! {
    /// Flags `Result<Result<T, E1>, E2>` in function signatures and type aliases.
    /// Almost always a mistake (`.map()` instead of `.and_then()`) or simplifiable.
    pub RESULT_RESULT,
    Warn,
    "nested `Result<Result<_, _>, _>` -- consider flattening into a single Result"
}

const MSG: &str =
    "nested `Result<Result<_, _>, _>` -- consider flattening into a single Result";
const HELP: &str =
    "use `.and_then()` to chain fallible operations, or unify the error types into a single enum";

/// Returns `true` if `ty` is `Result<Result<_, _>, _>` using `DefId` resolution.
///
/// Assumption: generic type params (e.g. `T` in `Result<T, E>`) are opaque
/// and won't match -- no false positives from `fn wrap<T>(v: T) -> Result<T, E>`.
fn is_nested_result<'tcx>(cx: &LateContext<'tcx>, ty: ty::Ty<'tcx>) -> bool {
    let ty::Adt(outer_adt, outer_args) = ty.kind() else {
        return false;
    };
    if !cx.tcx.is_diagnostic_item(sym::Result, outer_adt.did()) {
        return false;
    }

    let ok_ty = outer_args.type_at(0);
    let ty::Adt(inner_adt, _) = ok_ty.kind() else {
        return false;
    };
    cx.tcx.is_diagnostic_item(sym::Result, inner_adt.did())
}

pub struct ResultResult;

impl ResultResult {
    pub const fn new() -> Self {
        Self
    }
}

rustc_session::impl_lint_pass!(ResultResult => [RESULT_RESULT]);

impl<'tcx> LateLintPass<'tcx> for ResultResult {
    fn check_fn(
        &mut self,
        cx: &LateContext<'tcx>,
        kind: FnKind<'tcx>,
        decl: &'tcx FnDecl<'tcx>,
        _body: &'tcx Body<'tcx>,
        span: Span,
        def_id: rustc_hir::def_id::LocalDefId,
    ) {
        // Trait impl methods: signature is dictated by the trait, not the impl.
        if span.from_expansion()
            || matches!(kind, FnKind::Closure)
            || is_def_id_trait_method(cx, def_id)
        {
            return;
        }

        let ret_ty = return_ty(cx, rustc_hir::OwnerId { def_id });

        if is_nested_result(cx, ret_ty) {
            span_lint_and_help(cx, RESULT_RESULT, decl.output.span(), MSG, None, HELP);
        }
    }

    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        if item.span.from_expansion() {
            return;
        }

        let ItemKind::TyAlias(..) = &item.kind else {
            return;
        };

        let ty = cx.tcx.type_of(item.owner_id.def_id).instantiate_identity().skip_norm_wip();

        if is_nested_result(cx, ty) {
            span_lint_and_help(cx, RESULT_RESULT, item.span, MSG, None, HELP);
        }
    }
}
