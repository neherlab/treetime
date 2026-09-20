use clippy_utils::diagnostics::span_lint_and_then;
use clippy_utils::paths::{PathNS, lookup_path_str};
use clippy_utils::source::snippet_with_applicability;
use rustc_errors::Applicability;
use rustc_hir::def::Res;
use rustc_hir::def_id::DefId;
use rustc_hir::{Expr, ExprKind, GenericArg, QPath};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::def_id::LOCAL_CRATE;
use rustc_span::{Span, Symbol, sym};

rustc_session::declare_lint! {
    pub VALUE_CLONED_TO_DESERIALIZE,
    Warn,
    "`from_value` on a cloned `Value` -- deserialize from the borrowed value"
}

pub struct ValueClonedToDeserialize {
  from_value: Vec<DefId>,
  value: Vec<DefId>,
  serde_is_direct_dependency: bool,
}

impl ValueClonedToDeserialize {
  pub const fn new() -> Self {
    Self {
      from_value: Vec::new(),
      value: Vec::new(),
      serde_is_direct_dependency: false,
    }
  }
}

rustc_session::impl_lint_pass!(ValueClonedToDeserialize => [VALUE_CLONED_TO_DESERIALIZE]);

impl<'tcx> LateLintPass<'tcx> for ValueClonedToDeserialize {
  fn check_crate(&mut self, cx: &LateContext<'tcx>) {
    self.from_value = lookup_path_str(cx.tcx, PathNS::Value, "serde_json::from_value");
    self.value = lookup_path_str(cx.tcx, PathNS::Type, "serde_json::Value");
    let serde = Symbol::intern("serde");
    self.serde_is_direct_dependency = cx.tcx.crates(()).iter().any(|&krate| {
      cx.tcx.crate_name(krate) == serde
        && cx
          .tcx
          .extern_crate(krate)
          .is_some_and(|extern_crate| extern_crate.dependency_of == LOCAL_CRATE)
    });
  }

  fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
    if expr.span.from_expansion() || self.from_value.is_empty() {
      return;
    }
    let ExprKind::Call(callee, [arg]) = expr.kind else {
      return;
    };
    let ExprKind::Path(qpath) = &callee.kind else {
      return;
    };
    let Res::Def(_, did) = cx.qpath_res(qpath, callee.hir_id) else {
      return;
    };
    if !self.from_value.contains(&did) {
      return;
    }
    let ExprKind::MethodCall(method, receiver, [], _) = arg.kind else {
      return;
    };
    if method.ident.name != sym::clone {
      return;
    }
    let receiver_ty = cx.typeck_results().expr_ty(receiver);
    let ty::Adt(adt, _) = receiver_ty.peel_refs().kind() else {
      return;
    };
    if !self.value.contains(&adt.did()) {
      return;
    }
    let mut applicability = if self.serde_is_direct_dependency {
      Applicability::MachineApplicable
    } else {
      Applicability::MaybeIncorrect
    };
    let receiver_snippet = snippet_with_applicability(cx, receiver.span, "..", &mut applicability);
    let borrowed = if receiver_ty.is_ref() {
      receiver_snippet.into_owned()
    } else {
      format!("&{receiver_snippet}")
    };
    let replacement = match turbofish_type(qpath) {
      Some(target) => {
        let target = snippet_with_applicability(cx, target, "_", &mut applicability);
        format!("<{target} as serde::Deserialize>::deserialize({borrowed})")
      },
      None => format!("serde::Deserialize::deserialize({borrowed})"),
    };
    span_lint_and_then(
      cx,
      VALUE_CLONED_TO_DESERIALIZE,
      expr.span,
      "`from_value` copies the whole `Value` before parsing it",
      |diag| {
        diag.span_suggestion(
          expr.span,
          "deserialize from the borrowed value; `&Value` is a `Deserializer`",
          replacement,
          applicability,
        );
      },
    );
  }
}

fn turbofish_type(qpath: &QPath<'_>) -> Option<Span> {
  let QPath::Resolved(_, path) = qpath else {
    return None;
  };
  let args = path.segments.last()?.args?;
  let [GenericArg::Type(ty)] = args.args else {
    return None;
  };
  Some(ty.span)
}
