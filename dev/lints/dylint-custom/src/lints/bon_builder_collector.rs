use rustc_ast::{Item, ItemKind, MetaItem, MetaItemInner};
use rustc_lint::{EarlyContext, EarlyLintPass};
use rustc_span::symbol::sym;

use super::{DeriveInfo, STRUCT_DERIVES};

rustc_session::declare_lint! {
    /// Internal lint used to collect derive information from structs.
    ///
    /// Runs as a pre-expansion `EarlyLintPass` so it can inspect derive
    /// attributes before they are consumed by macro expansion.  Populates
    /// [`STRUCT_DERIVES`](super::STRUCT_DERIVES) with all derive trait names
    /// per struct (consumed by [`has_bon_builder`](super::has_bon_builder)
    /// and [`has_any_derive`](super::has_any_derive)).
    ///
    /// See [`has_bon_builder`](super::has_bon_builder) for the implications
    /// of name-only matching.
    pub BON_BUILDER_COLLECTOR,
    Allow,
    "internal: collects struct derive attributes"
}

/// Builds a [`DeriveInfo`] from the path of one derived trait, recording both
/// the last path segment and whether the path is exactly `bon::Builder`.
fn derive_info(meta: &MetaItem) -> Option<DeriveInfo> {
    let segments = &meta.path.segments;
    let last = segments.last()?;
    let last_segment = last.ident.name;
    // Path-aware `bon::Builder` detection: the second-to-last segment must be
    // `bon` so that `derive_builder::Builder` (same last segment) is excluded.
    let is_bon_builder = last_segment.as_str() == "Builder"
        && segments
            .iter()
            .rev()
            .nth(1)
            .is_some_and(|seg| seg.ident.name.as_str() == "bon");
    Some(DeriveInfo {
        last_segment,
        is_bon_builder,
    })
}

/// Records every trait listed inside a `derive(...)` meta-item list.
fn collect_from_derive_list(list: &[MetaItemInner], out: &mut Vec<DeriveInfo>) {
    for item in list {
        let MetaItemInner::MetaItem(meta) = item else {
            continue;
        };
        if let Some(info) = derive_info(meta) {
            out.push(info);
        }
    }
}

/// Returns `true` if a `cfg_attr` condition is unconditionally true -- the empty
/// `all()` predicate, which always holds. Conditional predicates (`test`,
/// `feature = "x"`, `not(..)`, etc.) are treated as not-always-applied, so the
/// derives they guard are not recorded.
fn is_always_true_cfg(cond: &MetaItemInner) -> bool {
    let MetaItemInner::MetaItem(meta) = cond else {
        return false;
    };
    meta.has_name(sym::all) && meta.meta_item_list().is_some_and(|list| list.is_empty())
}

/// Extracts derive info for every trait in `#[derive(...)]` *and*
/// `#[cfg_attr(<cond>, derive(...))]` attributes.
///
/// For `cfg_attr` the derives are only recorded when `<cond>` is
/// unconditionally true (the empty `all()` predicate). A conditional like
/// `cfg_attr(test, derive(bon::Builder))` or
/// `cfg_attr(feature = "x", derive(bon::Builder))` may not apply in this build,
/// so recording it would make downstream lints misfire on a struct that never
/// actually derives the builder.
fn collect_derive_names(attrs: &[rustc_ast::Attribute]) -> Vec<DeriveInfo> {
    let mut names = Vec::new();
    for attr in attrs {
        if attr.has_name(sym::derive) {
            if let Some(list) = attr.meta_item_list() {
                collect_from_derive_list(&list, &mut names);
            }
        } else if attr.has_name(sym::cfg_attr) {
            // `cfg_attr(<cond>, attr1, attr2, ...)`: the first element is the
            // condition; only honor the wrapped derives when it is always true.
            let Some(list) = attr.meta_item_list() else {
                continue;
            };
            let mut items = list.into_iter();
            let Some(cond) = items.next() else {
                continue;
            };
            if !is_always_true_cfg(&cond) {
                continue;
            }
            // Inspect every remaining inner attr, descending into any that are
            // themselves `derive(...)` lists.
            for item in items {
                let MetaItemInner::MetaItem(meta) = item else {
                    continue;
                };
                if meta.has_name(sym::derive)
                    && let Some(inner) = meta.meta_item_list()
                {
                    collect_from_derive_list(inner, &mut names);
                }
            }
        }
    }
    names
}

pub struct BonBuilderCollector;

rustc_session::impl_lint_pass!(BonBuilderCollector => [BON_BUILDER_COLLECTOR]);

impl EarlyLintPass for BonBuilderCollector {
    fn check_item(&mut self, _cx: &EarlyContext<'_>, item: &Item) {
        let ItemKind::Struct(ident, ..) = &item.kind else {
            return;
        };
        let derives = collect_derive_names(&item.attrs);
        if derives.is_empty() {
            return;
        }
        STRUCT_DERIVES.with(|map| {
            map.borrow_mut()
                .entry(ident.name)
                .or_default()
                .extend(derives);
        });
    }
}
