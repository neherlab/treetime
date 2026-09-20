use crate::adt_facts::{field_ty, has_fixed_repr, has_positional_fields};
use crate::baseline::emit;
use rustc_hir::{Item, ItemKind};
use rustc_lint::{LateContext, LateLintPass};

use crate::MordantConfig;

rustc_session::declare_lint! {
    /// Flags a named-field struct with `bool-cluster-min-bools` or more
    /// `bool` fields, however they are assigned. `n` bools are `2^n`
    /// possible combinations. If only some are valid, an enum names those.
    ///
    /// Silent on any struct with an explicit `repr` — its layout is dictated
    /// from outside Rust (a hardware register, a wire format), so all `2^n`
    /// states may genuinely be reachable.
    ///
    /// Runs only with `bool-cluster-enabled = true` in `dylint.toml`. Most
    /// structs it names are option bags whose states are all legal, so it is
    /// a survey to run once over a codebase, not a gate to keep on.
    pub BOOL_CLUSTER,
    Warn,
    "struct with several independent bool fields"
}

pub struct BoolCluster {
    pub config: &'static MordantConfig,
}

rustc_session::impl_lint_pass!(BoolCluster => [BOOL_CLUSTER]);

impl<'tcx> LateLintPass<'tcx> for BoolCluster {
    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        let ItemKind::Struct(..) = item.kind else {
            return;
        };
        let did = item.owner_id.to_def_id();
        let adt = cx.tcx.adt_def(did);
        if has_fixed_repr(adt) {
            return;
        }
        let variant = adt.non_enum_variant();
        if has_positional_fields(variant) {
            return;
        }
        let bools: Vec<_> = variant
            .fields
            .iter()
            .filter(|f| field_ty(cx, f).is_bool())
            .map(|f| f.name)
            .collect();
        if bools.len() < self.config.bool_cluster_min_bools {
            return;
        }
        let names: Vec<String> = bools.iter().map(|s| format!("`{s}`")).collect();
        let n = states(bools.len());
        let path = cx.tcx.def_path_str(did);
        emit(
            cx,
            BOOL_CLUSTER,
            cx.tcx.def_span(did),
            format!(
                "`{path}` has {} bool fields ({}). That is {n} possible combinations",
                bools.len(),
                names.join(", "),
            ),
            format!(
                "if only some combinations are valid, replace the bools with an enum of those. If the layout is fixed from outside, add a `repr` to `{path}` to say so"
            ),
        );
    }
}

fn states(n: usize) -> String {
    match u32::try_from(n).ok().and_then(|n| 2u64.checked_pow(n)) {
        Some(n) => n.to_string(),
        None => format!("2^{n}"),
    }
}

#[cfg(test)]
mod tests {
    use super::states;

    #[test]
    fn states_is_two_to_the_n() {
        assert_eq!(states(3), "8");
        assert_eq!(states(4), "16");
        assert_eq!(states(5), "32");
    }

    /// The count cannot be raised to a power of two in `u64`, so the message
    /// must say so rather than print a wrong number.
    #[test]
    fn states_does_not_overflow() {
        assert_eq!(states(64), "2^64");
        assert_eq!(states(usize::MAX), format!("2^{}", usize::MAX));
    }
}
