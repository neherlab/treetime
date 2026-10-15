# Branch-length validation tests match error substrings

The tests in `packages/treetime/src/optimize/__tests__/test_branch_length_validation.rs` check rejected branch lengths with `expect_err(..)` followed by `assert!(error.to_string().contains(..))` [packages/treetime/src/optimize/**tests**/test_branch_length_validation.rs#L46-L47](../../packages/treetime/src/optimize/__tests__/test_branch_length_validation.rs#L46-L47), [#L60-L61](../../packages/treetime/src/optimize/__tests__/test_branch_length_validation.rs#L60-L61), [#L69-L70](../../packages/treetime/src/optimize/__tests__/test_branch_length_validation.rs#L69-L70). A substring check passes when the message changes around the matched words, and a failure prints only `false` instead of the actual message.

## Proposed solution

Assert the exact error message with `assert_error!()` of `treetime_utils::testing`.
