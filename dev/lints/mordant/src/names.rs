//! The families the lints of this pack belong to. `GROUPS` drives
//! `LintStore::register_group` and `group:<name>` entries in `disabled`.

/// Every lint of the pack, in exactly one family. A family is a lint group
/// to rustc (`-A mordant_naming`, `#![warn(mordant_state)]`; see `group_id`)
/// and a `group:naming` entry in `disabled`.
pub const GROUPS: &[(&str, &[&str])] = &[
    (
        "state",
        &[
            "options_as_enum",
            "parallel_bools",
            "bool_cluster",
            "runtime_typestate",
            "always_unwrapped_option",
            "derived_field",
            "field_valid_only_when",
            "bool_beside_option",
            "parallel_vecs",
            "parallel_params",
            "stringly_state",
            "tuple_wants_struct",
            "some_still_unchecked",
        ],
    ),
    (
        "checks",
        &[
            "unchecked_construction",
            "defaulted_failure",
            "unchecked_input_len",
            "guard_blind_to_action",
            "stale_across_reentry",
            "error_collapsed_to_bool",
            "narrowed_two_ways",
            "cast_bypasses_from",
            "sentinel_integer",
        ],
    ),
    (
        "errors",
        &[
            "stringly_error",
            "stringified_error",
            "discarded_error",
            "unread_error_variant",
        ],
    ),
    (
        "enums",
        &[
            "wildcard_over_own_enum",
            "param_wider_than_callers",
            "return_wider_than_body",
        ],
    ),
    (
        "duplication",
        &[
            "same_match_twice",
            "reimplemented_helper",
            "generic_body_not_generic",
        ],
    ),
    (
        "naming",
        &[
            "bare_bool_args",
            "arg_named_like_other_param",
            "interchangeable_aliases",
            "index_of_other_kind",
            "unit_mismatch",
        ],
    ),
    (
        "keys_locks",
        &["key_not_identity", "insert_then_unwrap", "lock_order"],
    ),
    ("comments", &["stale_safety_comment", "stale_panic_message"]),
    ("unused", &["unused_pub"]),
    ("custom", &["forbidden_reach"]),
];

/// The group's name to rustc. Every library dylint loads shares one flat
/// namespace of lints and groups with rustc itself, and a bare `errors` or
/// `naming` in it would read as anyone's, so the id carries the pack's name.
pub fn group_id(group: &str) -> String {
    format!("mordant_{group}")
}

/// The lints a `disabled` entry spelled `group:<name>` stands for; `None`
/// for any other entry, including a `group:` naming no family.
pub fn group_members(entry: &str) -> Option<&'static [&'static str]> {
    let name = entry.strip_prefix("group:")?;
    GROUPS
        .iter()
        .find(|(g, _)| *g == name)
        .map(|(_, members)| *members)
}

#[cfg(test)]
mod tests {
    use super::{GROUPS, group_members};

    #[test]
    fn group_names_are_unique() {
        let names: Vec<&str> = GROUPS.iter().map(|(g, _)| *g).collect();
        for (i, g) in names.iter().enumerate() {
            assert!(!names[i + 1..].contains(g), "{g} twice");
        }
    }

    #[test]
    fn group_members_reads_only_the_group_spelling() {
        assert_eq!(
            group_members("group:duplication"),
            Some(
                &[
                    "same_match_twice",
                    "reimplemented_helper",
                    "generic_body_not_generic",
                ][..]
            )
        );
        assert_eq!(group_members("group:nope"), None);
        assert_eq!(group_members("duplication"), None);
        assert_eq!(group_members("lock_order"), None);
    }
}
