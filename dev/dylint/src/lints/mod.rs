pub mod bon_builder_collector;
pub mod debug_remnants;
pub mod fallible_new;
mod hir_refs;
pub mod needless_builder;
pub mod panic_in_drop;
pub mod prefer_error_macros;
pub mod proper_error_type;
pub mod result_result;
pub mod suggest_builder;
mod suppression;
pub mod topological_ordering;
pub mod unclear_exports;

use core::cell::RefCell;
use std::collections::{HashMap, HashSet};

use rustc_span::Symbol;

/// A single `#[derive(...)]` entry recorded during the pre-expansion pass.
///
/// We keep both the *last* path segment (for name-only skip-derive matching,
/// e.g. `Default` matching both `Default` and `std::default::Default`) and a
/// path-aware flag identifying `bon::Builder` specifically, so that
/// `derive_builder::Builder`, which shares the last segment `Builder`, is not
/// confused for it.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct DeriveInfo {
    /// The last path segment of the derive (e.g. `Builder`, `Default`).
    pub last_segment: Symbol,
    /// Whether the derive path is exactly `bon::Builder`.
    pub is_bon_builder: bool,
}

thread_local! {
    /// Maps struct names to the set of derives found on them during the
    /// pre-expansion pass.  Populated by [`BonBuilderCollector`] and consumed
    /// via [`has_any_derive`] / [`has_bon_builder`].
    pub static STRUCT_DERIVES: RefCell<HashMap<Symbol, HashSet<DeriveInfo>>> = RefCell::new(HashMap::new());
}

/// Returns `true` if any of the given derive names were found (by *last path
/// segment*) on a struct with the given name during the pre-expansion pass.
///
/// Matching is by last segment so that, e.g., `Default` matches both
/// `#[derive(Default)]` and `#[derive(std::default::Default)]`.
pub fn has_any_derive(name: Symbol, derives: &[Symbol]) -> bool {
    STRUCT_DERIVES.with(|map| {
        map.borrow().get(&name).is_some_and(|set| {
            set.iter()
                .any(|info| derives.contains(&info.last_segment))
        })
    })
}

/// Returns `true` if a struct with the given name was found to have
/// `#[derive(bon::Builder)]` (path-aware) during the pre-expansion pass.
///
/// Matching is path-aware: `derive_builder::Builder` shares the last segment
/// `Builder` but is *not* `bon::Builder`, so it is excluded.
///
/// **Limitation:** Uses name-only *keying* (not `DefId`) because the
/// pre-expansion AST pass runs before name resolution.  If two structs in
/// different modules share the same name and only one derives `bon::Builder`,
/// both will be treated identically.  Switching to a `LateLintPass` would fix
/// this at the cost of not seeing derives consumed by macro expansion.
pub fn has_bon_builder(name: Symbol) -> bool {
    STRUCT_DERIVES.with(|map| {
        map.borrow()
            .get(&name)
            .is_some_and(|set| set.iter().any(|info| info.is_bon_builder))
    })
}
