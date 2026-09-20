//! The rev the loaded pack was built from (T-328), surfaced through the
//! same channel every real lint here already uses to describe itself.
//!
//! `cargo dylint list` does not run rustc's own `-W help`: the driver
//! snapshots `LintStore::get_lints()` before and after each library's
//! `register_lints` runs, prints `name  level  desc` for every lint that
//! appeared in the diff, then exits (`dylint_driver::list_lints`). A
//! `declare_lint!` is a `pub static … : &Lint`, so registering one with the
//! rev baked into its `desc` makes that line dylint's own answer to "what
//! rev is this" — no separate query, no parsing a build artefact.
//!
//! This is not a check and never fires: no pass reads `MORDANT_SOURCE_REV`,
//! so it is registered directly here rather than through `register()`'s
//! `Registrar`, and deliberately outside `names::GROUPS` — a group implies
//! an invariant a caller might `-A`/`disabled` away, and there is none to
//! turn off. Keeping it out of `Registrar::known` also keeps it out of
//! `every_lint_is_in_exactly_one_group_and_the_group_resolves_to_it`, whose
//! whole point is to keep the *checking* lints' bookkeeping exhaustive.

rustc_session::declare_lint! {
    /// Never fires. Registered only so its `desc` — the git rev this build
    /// was compiled from — appears in `cargo dylint list`'s output.
    pub MORDANT_SOURCE_REV,
    Allow,
    concat!("mordant: pack built from rev ", env!("MORDANT_SOURCE_REV"))
}

/// Registers the marker on the raw store. Called from `register_lints`
/// (the FFI entry point), not from `register()`, so `registered_store()` in
/// `lib.rs`'s tests — which exercises `register()` alone — never sees it.
pub fn register(s: &mut rustc_lint::LintStore) {
    s.register_lints(&[MORDANT_SOURCE_REV]);
}

#[cfg(test)]
mod tests {
    use super::register;

    /// `cargo dylint list` reads exactly two fields off a registered lint:
    /// its name and its `desc`. This is the contract those have to meet —
    /// findable under the name, and carrying a rev that is not empty.
    #[test]
    fn register_adds_a_lint_whose_desc_carries_a_nonempty_rev() {
        let mut store = rustc_lint::LintStore::new();
        register(&mut store);
        let lint = store
            .get_lints()
            .iter()
            .find(|l| l.name_lower() == "mordant_source_rev")
            .expect("mordant_source_rev registered");
        let rev = lint
            .desc
            .rsplit(' ')
            .next()
            .expect("desc has at least one word");
        assert!(!rev.is_empty(), "desc {:?} carries no rev", lint.desc);
        assert_ne!(
            rev, "rev",
            "desc {:?} looks unstamped (trailing word is the label, not a value)",
            lint.desc
        );
    }
}
