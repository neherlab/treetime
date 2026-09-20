#![feature(rustc_private)]
#![feature(default_field_values)]
#![warn(unused_extern_crates)]

extern crate rustc_abi;
extern crate rustc_ast;
extern crate rustc_data_structures;
extern crate rustc_errors;
extern crate rustc_hir;
extern crate rustc_index;
extern crate rustc_lint;
extern crate rustc_metadata;
extern crate rustc_middle;
extern crate rustc_mir_dataflow;
extern crate rustc_session;
extern crate rustc_span;

dylint_linting::dylint_library!();

use rustc_data_structures::sync;

mod adt_facts;
mod always_unwrapped_option;
mod arg_named_like_other_param;
mod bare_bool_args;
mod baseline;
mod bool_beside_option;
mod bool_cluster;
mod cast_bypasses_from;
mod claims;
mod ctor_flow;
mod defaulted_failure;
mod derived_field;
mod discarded_error;
mod enum_facts;
mod error_collapsed_to_bool;
mod field_valid_only_when;
mod forbidden_reach;
mod generic_body_not_generic;
mod guard_blind_to_action;
mod hir_clone;
mod hir_shapes;
mod index_of_other_kind;
mod insert_then_unwrap;
mod interchangeable_aliases;
mod key_not_identity;
mod lock_order;
mod mir_flow;
mod names;
mod narrowed_two_ways;
mod options_as_enum;
mod parallel_bools;
mod parallel_params;
mod parallel_vecs;
mod param_wider_than_callers;
mod reimplemented_helper;
mod return_wider_than_body;
mod runtime_typestate;
mod same_match_twice;
mod sentinel_integer;
mod some_still_unchecked;
mod source_rev;
mod stale_across_reentry;
mod stale_panic_message;
mod stale_safety_comment;
mod stringified_error;
mod stringly_error;
mod stringly_state;
mod tuple_wants_struct;
mod unchecked_construction;
mod unchecked_input_len;
mod unit_mismatch;
mod unread_error_variant;
mod unused_pub;
mod variant_flow;
mod wildcard_over_own_enum;

/// Read from `dylint.toml` under `[mordant]` in the linted workspace root.
///
/// `bool_cluster` is allowed on this one struct, and it is the lawful-lattice
/// case the lint's own help describes rather than an exemption from it. The
/// bools are opt-ins belonging to *different* lints: every combination is
/// reachable from a `dylint.toml` and each means what it says, so there is no
/// invariant between them for a type to carry. The field set is also this
/// pack's public interface — every field is a TOML key, and Scarlet's `xtask`
/// gate reads these field names out of the pinned source to decide which keys
/// it may set — so grouping them into sub-structs would rename user-visible
/// keys to satisfy a lint about internal invariants.
///
/// `dylint_linting::config_or_default` returns `Default` when the linted
/// workspace has no `dylint.toml`; the container-level `serde(default)`
/// fills any omitted key from it too.
#[derive(Default, serde::Deserialize)]
#[cfg_attr(test, derive(Debug, PartialEq))]
#[serde(rename_all = "kebab-case", default, deny_unknown_fields)]
#[cfg_attr(dylint_lib = "mordant", allow(bool_cluster))]
pub struct MordantConfig {
    /// Lints, by name, that stay registered (so an `allow` of one still
    /// resolves) but never run. `group:<name>` stands for every lint in that
    /// family (`names::GROUPS`).
    pub disabled: Vec<String>,
    /// Fully qualified paths of types that are never a valid map key in this
    /// project (e.g. a span type with no file identity). Empty means silent.
    pub key_not_identity_types: Vec<String>,
    /// Opt-in key-expression forms: "to-bits", "ptr-cast". Both are legitimate
    /// in interning code, so neither is on by default.
    pub key_not_identity_forms: Vec<String>,
    /// Also flag `Box<dyn Error>` as a stringly error type.
    pub stringly_error_include_box_dyn: bool,
    /// Fully qualified method paths that never produce a valid map key (e.g. a
    /// NaN-boxing `Value::to_bits`). Checked in key position of `insert`/`entry`.
    pub key_not_identity_methods: Vec<String>,
    /// Opt-in: also flag composite keys (tuples, structs one level deep) that
    /// carry a denied type without one of the fixing types beside it.
    pub key_not_identity_composite: bool,
    /// Types whose presence in a composite key restores identity (e.g. the
    /// file id that gives a span a coordinate space).
    pub key_not_identity_fixes: Vec<String>,
    /// Error types that mean "the environment refused" (allocation, IO,
    /// syscall), added to the built-in std list. A constructor exit failing
    /// with one of these is never treated as validating a field.
    pub validator_resource_errors: Vec<String>,
    /// Minimum Option fields for `options_as_enum` to consider a struct.
    pub options_as_enum_min_fields: usize = 2,
    /// `wildcard_over_own_enum` stays silent above this many variants.
    pub wildcard_over_own_enum_max_variants: usize = 12,
    /// Bool fields at which `bool_cluster` fires.
    pub bool_cluster_min_bools: usize = 3,
    /// Construction sites at which `derived_field` will read a
    /// correspondence between two fields.
    pub derived_field_min_sites: usize = 2,
    /// Expression nodes below which `reimplemented_helper` does not compare
    /// a body, so one-line accessors and constructors never pair up.
    pub reimplemented_helper_min_nodes: usize = 12,
    /// Ratchet file name, resolved upward from each crate's manifest dir. Runs
    /// suppress up to the recorded count per (lint, file) and surface only new
    /// findings. Regenerate with `MORDANT_BASELINE_WRITE=1`.
    pub baseline: Option<String>,
    /// Reachability bans: from each matching root, no call path may reach a
    /// banned definition. Findings carry the witness path.
    pub forbidden_reach: Vec<forbidden_reach::ReachRule>,
    /// This project's own re-entry points, for `stale_across_reentry`, on top
    /// of the built-in closure / fn-pointer / `dyn` / `.await` set: paths
    /// matched by `::`-segment suffix (`Vm::run_callback`, `run_callback`;
    /// a method of a trait impl matches under its type or its trait,
    /// `Worker::run_job` or `Runner::run_job`), with a trailing `*` on the
    /// last segment matching by prefix (`dispatch*`).
    pub stale_across_reentry_callees: Vec<String>,
    /// Callees `defaulted_failure` takes as rejecting their argument without
    /// reading their body: parsers in other crates, or local ones whose
    /// failure is built by combinators. Spelled like
    /// `validator-resource-errors` (a full path, `crate::name`, or a bare
    /// name). Empty by default; the lint then reports only callees whose
    /// body it can see the check in.
    pub defaulted_failure_callees: Vec<String>,
    /// Opt-in: run `bool_cluster`. Off by default because most structs it
    /// names are option bags; the state machines among them are found by
    /// running it once, not by gating on it.
    pub bool_cluster_enabled: bool,
    /// Opt-in: run `stale_safety_comment`. Off by default because a name a
    /// crate cannot see is usually defined in C++, a script, or a downstream
    /// crate; the stale ones are found by running it once.
    pub stale_safety_comment_enabled: bool,
    /// Error types whose failure `defaulted_failure` does not count, on top
    /// of `validator-resource-errors`: errors that are already recorded
    /// somewhere else by the time they are returned (a "JS exception is
    /// pending" marker, a diagnostic already pushed to a log), so defaulting
    /// them hides nothing. Same spellings as `validator-resource-errors`.
    pub defaulted_failure_ignored_errors: Vec<String>,
    /// Opt-in: run `unchecked_input_len`. Off by default because most of
    /// what it names on a codebase whose callers vouch for their lengths is a
    /// value the function also uses as some other value's limit (TRIAGE.md),
    /// which nothing inside the function tells from a missed check; run it
    /// once over parsing code and read the list.
    pub unchecked_input_len_enabled: bool,
    /// Opt-in: run `parallel_params`. Off by default because a buffer and a
    /// cursor into it, or a precedence level and the flags in force at it,
    /// pass between functions together by design, and nothing in the
    /// signatures tells those from a value nobody declared; run it once and
    /// read the list.
    pub parallel_params_enabled: bool,
    /// Opt-in: run `some_still_unchecked`. Off by default because a `Some` that fails
    /// the check usually is meant to read as absent; the lint is a sweep for
    /// the places where a `.filter(..)` or a narrower type says so instead.
    pub some_still_unchecked_enabled: bool,
    /// Functions a parameter group must pass between, unchanged, before
    /// `parallel_params` names it.
    pub parallel_params_min_fns: usize = 3,
    /// Opt-in: run `generic_body_not_generic`. Off by default because its
    /// counts are exact but the bytes a fix saves are not. They depend on
    /// whether the compiler inlines the new function back and on whether the
    /// linker already folds identical copies. For small parts the saving
    /// measured near zero. It is a survey to run once over a size-sensitive
    /// crate, largest findings first.
    pub generic_body_not_generic_enabled: bool,
    /// The smallest number of statements the shared part of a generic fn's
    /// body needs before `generic_body_not_generic` names the fn. The shared
    /// part is one set of the body's blocks with a single entry and a single
    /// exit that mentions no type or const parameter, and whose values in and
    /// out have types free of those parameters. The count is of MIR statements
    /// and terminators; storage markers and plain jumps are not counted, so
    /// a fn that only forwards its arguments stays under it.
    ///
    /// There is no second threshold on how much of the body the shared part
    /// covers. A part that passes moves into a separate non-generic fn the
    /// same way however large the rest of the body is, and saves the same
    /// bytes: its own, and only if the compiler keeps that fn out of line.
    /// A threshold on its share of the body would only stop the lint from
    /// reporting fns whose fix is the same.
    pub generic_body_not_generic_min_statements: usize = 24,
    /// Distinct concrete generic-argument sets this crate must instantiate
    /// the fn at before its body counts as duplicated. One instantiation
    /// duplicates nothing.
    pub generic_body_not_generic_min_instantiations: usize = 2,
}

#[expect(clippy::no_mangle_with_rust_abi)]
#[unsafe(no_mangle)]
pub fn register_lints(sess: &rustc_session::Session, s: &mut rustc_lint::LintStore) {
    dylint_linting::init_config(sess);
    source_rev::register(s);
    let config: MordantConfig = dylint_linting::config_or_default(env!("CARGO_PKG_NAME"));
    let config: &'static MordantConfig = Box::leak(Box::new(config));
    baseline::setup(&config.baseline);
    let unknown = register(config, s);
    if !unknown.is_empty() {
        sess.dcx().warn(format!(
            "mordant: `disabled` in dylint.toml names no lint of this pack: {}",
            unknown.join(", ")
        ));
    }
}

/// Everything `register_lints` does to the store, apart from the session so
/// a test can run it against a bare `LintStore`. Returns the `disabled`
/// entries that named nothing.
fn register(config: &'static MordantConfig, s: &mut rustc_lint::LintStore) -> Vec<String> {
    use {
        always_unwrapped_option::AlwaysUnwrappedOption,
        arg_named_like_other_param::ArgNamedLikeOtherParam, bare_bool_args::BareBoolArgs,
        baseline::BaselineWriter, bool_cluster::BoolCluster, cast_bypasses_from::CastBypassesFrom,
        defaulted_failure::DefaultedFailure, derived_field::DerivedField,
        discarded_error::DiscardedError, error_collapsed_to_bool::ErrorCollapsedToBool,
        field_valid_only_when::FieldValidOnlyWhen, forbidden_reach::ForbiddenReach,
        generic_body_not_generic::GenericBodyNotGeneric, guard_blind_to_action::GuardBlindToAction,
        index_of_other_kind::IndexOfOtherKind, insert_then_unwrap::InsertThenUnwrap,
        interchangeable_aliases::InterchangeableAliases, key_not_identity::KeyNotIdentity,
        lock_order::LockOrder, narrowed_two_ways::NarrowedTwoWays, options_as_enum::OptionsAsEnum,
        parallel_bools::ParallelBools, parallel_params::ParallelParams,
        parallel_vecs::ParallelVecs, param_wider_than_callers::ParamWiderThanCallers,
        return_wider_than_body::ReturnWiderThanBody, runtime_typestate::RuntimeTypestate,
        sentinel_integer::SentinelInteger, stale_across_reentry::StaleAcrossReentry,
        stale_panic_message::StalePanicMessage, stale_safety_comment::StaleSafetyComment,
        stringified_error::StringifiedError, stringly_error::StringlyError,
        stringly_state::StringlyState, tuple_wants_struct::TupleWantsStruct,
        unchecked_construction::UncheckedConstruction, unchecked_input_len::UncheckedInputLen,
        unit_mismatch::UnitMismatch, unread_error_variant::UnreadErrorVariant,
        wildcard_over_own_enum::WildcardOverOwnEnum,
    };
    let disabled = resolve_disabled(&config.disabled);
    let mut r = Registrar {
        store: s,
        disabled: &disabled,
        known: Vec::new(),
    };
    r.add(true, move || StringlyError { config });
    r.add(true, move || KeyNotIdentity::new(config));
    r.add(true, || StringifiedError);
    r.add(true, move || OptionsAsEnum::new(config));
    r.add(true, ParallelBools::default);
    r.add(config.bool_cluster_enabled, move || BoolCluster { config });
    r.add(true, move || UncheckedConstruction::new(config));
    r.add(true, UnreadErrorVariant::default);
    r.add(true, GuardBlindToAction::default);
    r.add(
        config.stale_safety_comment_enabled,
        StaleSafetyComment::default,
    );
    r.add(true, || UnitMismatch);
    r.add(true, StalePanicMessage::default);
    r.add(true, LockOrder::default);
    r.add(true, move || ForbiddenReach::new(config));
    r.add(true, RuntimeTypestate::default);
    r.add(true, AlwaysUnwrappedOption::default);
    r.add(true, || InsertThenUnwrap);
    r.add(true, ParamWiderThanCallers::default);
    r.add(true, ReturnWiderThanBody::default);
    r.add(true, move || WildcardOverOwnEnum { config });
    r.add(true, || DiscardedError);
    r.add(true, move || DerivedField::new(config));
    r.add(true, move || StaleAcrossReentry { config });
    r.add(true, move || DefaultedFailure::new(config));
    r.add(config.unchecked_input_len_enabled, || UncheckedInputLen);
    r.add(true, || ArgNamedLikeOtherParam);
    r.add(true, move || CastBypassesFrom::new(config));
    r.add(true, same_match_twice::SameMatchTwice::default);
    r.add(true, move || {
        reimplemented_helper::ReimplementedHelper::new(config)
    });
    r.add(true, FieldValidOnlyWhen::default);
    r.add(true, ErrorCollapsedToBool::default);
    r.add(true, NarrowedTwoWays::default);
    r.add(true, || IndexOfOtherKind);
    r.add(true, ParallelVecs::default);
    r.add(true, bool_beside_option::BoolBesideOption::default);
    r.add(true, SentinelInteger::default);
    r.add(true, StringlyState::default);
    r.add(config.parallel_params_enabled, move || {
        ParallelParams::new(config)
    });
    r.add(true, BareBoolArgs::default);
    r.add(true, TupleWantsStruct::default);
    r.add(true, || InterchangeableAliases);
    r.add(config.some_still_unchecked_enabled, || {
        some_still_unchecked::SomeStillUnchecked
    });
    r.add(config.generic_body_not_generic_enabled, move || {
        GenericBodyNotGeneric::new(config)
    });
    // Last, so its check_crate_post flushes after every lint has recorded.
    r.add(true, unused_pub::UnusedPub::default);
    r.add(true, || BaselineWriter);
    r.groups(names::GROUPS);
    unknown_names(&disabled, &r.known)
}

/// Registers lints and passes through one seam, so a pass that is off (an
/// opt-in whose key is absent, or a lint listed under `disabled`) still has
/// its lint registered and `allow(..)` / `-A` of it still resolves.
struct Registrar<'a> {
    store: &'a mut rustc_lint::LintStore,
    disabled: &'a [String],
    /// Every lint name registered so far, to tell which `disabled` entries
    /// name nothing.
    known: Vec<String>,
}

impl Registrar<'_> {
    fn add<T: for<'tcx> rustc_lint::LateLintPass<'tcx> + 'static>(
        &mut self,
        enabled: bool,
        make: impl Fn() -> T + sync::DynSend + sync::DynSync + 'static,
    ) {
        let lints = make().get_lints();
        self.store.register_lints(&lints);
        let names: Vec<String> = lints.iter().map(|l| l.name_lower()).collect();
        let run = enabled && !all_disabled(&names, self.disabled);
        self.known.extend(names);
        if run {
            self.store.register_late_pass(move |_| Box::new(make()));
        }
    }

    /// After every `add`: each family becomes a lint group, so one `-A` /
    /// `#![allow(..)]` of `names::group_id` covers its members.
    fn groups(&mut self, table: &[(&str, &[&str])]) {
        for (group, members) in table {
            let ids = members
                .iter()
                .flat_map(|m| match self.store.find_lints(m) {
                    Some(ids) => ids.to_vec(),
                    None => panic!("mordant: group {group} lists {m}, which is not a lint"),
                })
                .collect();
            let id: &'static str = Box::leak(names::group_id(group).into_boxed_str());
            self.store.register_group(true, id, None, ids);
        }
    }
}

/// `disabled` as written in `dylint.toml`, with `group:<name>` expanded to
/// every lint in that family. An entry that names nothing is kept as written
/// for `unknown_names` to report.
fn resolve_disabled(written: &[String]) -> Vec<String> {
    let mut disabled = Vec::new();
    for entry in written {
        match names::group_members(entry) {
            Some(members) => disabled.extend(members.iter().map(|m| m.to_string())),
            None => disabled.push(entry.clone()),
        }
    }
    disabled
}

/// A pass is skipped when every lint it declares is disabled. A pass with
/// no lints (the baseline writer) always runs.
fn all_disabled(lints: &[String], disabled: &[String]) -> bool {
    !lints.is_empty() && lints.iter().all(|l| disabled.contains(l))
}

fn unknown_names(disabled: &[String], known: &[String]) -> Vec<String> {
    disabled
        .iter()
        .filter(|d| !known.contains(d))
        .cloned()
        .collect()
}

#[test]
fn ui() {
    dylint_testing::ui::Test::src_base(env!("CARGO_PKG_NAME"), "ui")
        .dylint_toml(
            r#"
            [mordant]
            # Its fixtures are full of `pub` items nothing calls; it has its own suite.
            disabled = ["unused_pub"]
            key-not-identity-types = ["Span"]
            key-not-identity-forms = ["to-bits", "ptr-cast"]
            key-not-identity-methods = ["Value::to_raw"]
            key-not-identity-composite = true
            key-not-identity-fixes = ["FileId"]
            stale-across-reentry-callees = ["Vm::run_callback", "dispatch*", "Worker::run_job", "Runner::schedule"]
            defaulted-failure-callees = ["from_str_radix", "listed_by_config"]
            defaulted-failure-ignored-errors = ["Pending"]
            bool-cluster-enabled = true
            stale-safety-comment-enabled = true
            unchecked-input-len-enabled = true
            parallel-params-enabled = true
            some-still-unchecked-enabled = true
            generic-body-not-generic-enabled = true

            [[mordant.forbidden-reach]]
            from = "hot_path"
            never = ["std::vec::Vec::push"]

            [[mordant.forbidden-reach]]
            from = "two_bans"
            never = ["std::vec::Vec::push", "Option::expect"]

            [[mordant.forbidden-reach]]
            from = "one_ban_twice"
            never = ["std::vec::Vec::push"]

            [[mordant.forbidden-reach]]
            from = "index_root"
            never = ["panic_bounds_check"]

            [[mordant.forbidden-reach]]
            from = "add_overflow_root"
            never = ["panic_const_add_overflow"]

            [[mordant.forbidden-reach]]
            from = "div_zero_root"
            never = ["panic_const_div_by_zero"]

            [[mordant.forbidden-reach]]
            from = "rem_zero_root"
            never = ["panic_const_rem_by_zero"]

            [[mordant.forbidden-reach]]
            from = "neg_overflow_root"
            never = ["panic_const_neg_overflow"]

            # The two controls: a live ban on a family their bodies do not
            # reach, so a silent run here is the lint discriminating and not
            # the rule failing to resolve.
            [[mordant.forbidden-reach]]
            from = "wrong_family_root"
            never = ["panic_const_add_overflow"]

            [[mordant.forbidden-reach]]
            from = "no_assert_root"
            never = ["panic_const_add_overflow"]
            "#,
        )
        .run();
}

/// The `ui` fixtures for the opt-in lints run with their keys on; these
/// re-run the same shapes with the keys absent and expect nothing.
#[test]
fn ui_opt_in_lints_are_off_without_their_key() {
    dylint_testing::ui::Test::src_base(env!("CARGO_PKG_NAME"), "ui_off")
        .dylint_toml("[mordant]\ndisabled = [\"unused_pub\"]\n")
        .run();
}

/// `unused_pub` alone, since every other fixture is made of unused items.
#[test]
fn ui_unused_pub() {
    dylint_testing::ui::Test::src_base(env!("CARGO_PKG_NAME"), "ui_unused_pub")
        .dylint_toml("[mordant]\n")
        .run();
}

/// `config_or_default` returns `Default` when the linted workspace has no
/// `dylint.toml`. A threshold that lost its `= N` would default to 0, which
/// turns `wildcard_over_own_enum` off (`n > 0` for every enum) and makes
/// `options_as_enum` consider every struct.
#[test]
fn config_default_thresholds_match_docs() {
    let c = MordantConfig::default();
    assert_eq!(c.options_as_enum_min_fields, 2);
    assert_eq!(c.wildcard_over_own_enum_max_variants, 12);
    assert_eq!(c.bool_cluster_min_bools, 3);
    assert_eq!(c.derived_field_min_sites, 2);
    assert_eq!(c.reimplemented_helper_min_nodes, 12);
    assert!(!c.bool_cluster_enabled);
    assert!(!c.stale_safety_comment_enabled);
    assert!(!c.unchecked_input_len_enabled);
    assert!(!c.parallel_params_enabled);
    assert!(!c.some_still_unchecked_enabled);
    assert!(!c.generic_body_not_generic_enabled);
    assert_eq!(c.parallel_params_min_fns, 3);
    assert_eq!(c.generic_body_not_generic_min_statements, 24);
    assert_eq!(c.generic_body_not_generic_min_instantiations, 2);
}

/// An empty table (file present, keys omitted) must not drift from
/// `Default`.
#[test]
fn config_omitted_toml_keys_use_the_same_thresholds() {
    let parsed: MordantConfig = toml::from_str("").expect("empty document is an empty table");
    assert_eq!(parsed, MordantConfig::default());
}

#[test]
fn config_explicit_zero_thresholds_are_honored() {
    let parsed: MordantConfig = toml::from_str(
        "options-as-enum-min-fields = 0\n\
         wildcard-over-own-enum-max-variants = 0\n\
         bool-cluster-min-bools = 0\n",
    )
    .expect("explicit zeros parse");
    assert_eq!(parsed.options_as_enum_min_fields, 0);
    assert_eq!(parsed.wildcard_over_own_enum_max_variants, 0);
    assert_eq!(parsed.bool_cluster_min_bools, 0);
}

/// A misspelled key must not vanish into `default`: `key_not_identity_types`'s real spelling
/// is `key-not-identity-types`, and two consumers have carried a different, wrong spelling of
/// two different keys in production with no error (`language`'s `nonidentity-key-types` and
/// `scarlet`'s `wildcard-local-enum-max-variants`) — this is the shape that let both happen.
#[test]
fn unknown_key_is_an_error_naming_it() {
    let err = toml::from_str::<MordantConfig>("wildcard-local-enum-max-variants = 64\n")
        .expect_err("a key matching no field must not parse");
    let msg = err.to_string();
    assert!(
        msg.contains("wildcard-local-enum-max-variants"),
        "error should name the offending key: {msg}"
    );
}

#[test]
fn config_disabled_parses_and_defaults_empty() {
    assert!(MordantConfig::default().disabled.is_empty());
    let parsed: MordantConfig =
        toml::from_str("disabled = [\"runtime_typestate\", \"lock_order\"]\n")
            .expect("disabled parses");
    assert_eq!(parsed.disabled, ["runtime_typestate", "lock_order"]);
}

#[test]
fn disabled_skips_a_pass_only_when_every_lint_it_declares_is_named() {
    let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    assert!(all_disabled(
        &names(&["runtime_typestate"]),
        &names(&["runtime_typestate"])
    ));
    assert!(!all_disabled(
        &names(&["runtime_typestate"]),
        &names(&["lock_order"])
    ));
    assert!(!all_disabled(&names(&["a", "b"]), &names(&["a"])));
    assert!(!all_disabled(&[], &names(&["runtime_typestate"])));
}

#[test]
fn disabled_names_that_match_no_lint_are_reported() {
    let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let disabled = names(&[
        "runtime_typestate",
        "runtime_typstate",
        "clippy::unwrap_used",
    ]);
    let known = names(&["runtime_typestate", "lock_order"]);
    assert_eq!(
        unknown_names(&disabled, &known),
        ["runtime_typstate", "clippy::unwrap_used"]
    );
}

#[cfg(test)]
fn registered_store(config: MordantConfig) -> (rustc_lint::LintStore, Vec<String>) {
    let mut store = rustc_lint::LintStore::new();
    let unknown = register(Box::leak(Box::new(config)), &mut store);
    (store, unknown)
}

#[test]
fn every_lint_is_in_exactly_one_group_and_the_group_resolves_to_it() {
    let (store, _) = registered_store(MordantConfig::default());
    for lint in store.get_lints() {
        let name = lint.name_lower();
        let homes: Vec<&str> = names::GROUPS
            .iter()
            .filter(|(_, members)| members.contains(&name.as_str()))
            .map(|(g, _)| *g)
            .collect();
        assert_eq!(homes.len(), 1, "{name} is in {homes:?}");
        let group = store
            .find_lints(&names::group_id(homes[0]))
            .expect("group is registered");
        assert!(group.contains(&rustc_lint::LintId::of(lint)), "{name}");
    }
    let grouped: usize = names::GROUPS.iter().map(|(_, m)| m.len()).sum();
    assert_eq!(grouped, store.get_lints().len());
}

#[test]
fn disabled_expands_a_group_to_its_members() {
    let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let disabled = resolve_disabled(&names(&["group:duplication", "group:nope"]));
    assert_eq!(
        disabled,
        [
            "same_match_twice",
            "reimplemented_helper",
            "generic_body_not_generic",
            "group:nope"
        ]
    );
    let (_, unknown) = registered_store(MordantConfig {
        disabled: names(&["group:duplication", "group:nope"]),
        ..MordantConfig::default()
    });
    assert_eq!(unknown, ["group:nope"]);
}
