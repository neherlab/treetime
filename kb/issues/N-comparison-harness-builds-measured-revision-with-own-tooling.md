# Comparison harness builds the measured revision with its own tooling

> [!IMPORTANT]
> **Decision required.** The containment owner of the baseline build is not approved. A trusted harness revision (O1) keeps containment outside the measured revision; the baseline's own build tooling (O2) matches how that revision builds but lets it control containment. See "Decision axes".

The output-equality harness `dev/smoke` checks out the baseline (by default `rust`) into a worktree and builds it by running that revision's `./dev/docker/run` and its build task (`just br`, or `./dev/dev br` in older revisions). A revision under measurement therefore defines the mechanism intended to contain its own build execution.

`fn build_baseline_binary()` in `dev/smoke` adds a detached worktree at the baseline commit, and `fn build_in_checkout()` runs that checkout's own `./dev/docker/run` and build task [dev/smoke#L888-L924](../../dev/smoke#L888-L924). `WRAPPER` defaults to the relative path `./dev/docker/run`, and `run_checked()` runs it with `cwd=checkout`, so the baseline worktree's own wrapper runs. This makes the measured revision part of the build and containment mechanism. It also means that build settings `dev/smoke` passes to the baseline, such as `TREETIME_CARGO_HOME`, take effect only when the baseline revision's `dev/docker/run` reads them.

## Decision axes

### Containment owner

- O1. Use a separately identified trusted harness revision to build the baseline binary.
- O2. Execute the baseline revision's own build tooling. This matches how that revision builds but grants the measured revision control of containment.

Recommendation: O1. The baseline commit's build behavior can be inspected separately without executing it.

## Recommendation

Build the baseline binary from a trusted harness revision rather than the checked-out baseline's own scripts. Implementation waits until the containment-owner policy is approved.

## Proposed implementation plan

The plan stops at the open decision; implement one branch after approval.

- O1, trusted harness revision:
  1. In `dev/smoke`, resolve the wrapper to an absolute path in the current checkout (`PROJECT_ROOT / "dev/docker/run"`) when `SMOKE_WRAPPER` is unset, and run it with the baseline worktree as the mounted project directory. Check first how `dev/docker/run` picks the mounted directory (script location or working directory) and add an explicit argument or environment variable if it uses the script location
  2. Keep `release_build_cmd()` choosing the task runner of the baseline (`just br` or `./dev/dev br`), because the build recipe is part of what the revision builds; only the container boundary moves to the harness
  3. Delete the docstring sentence of `build_in_checkout()` about revisions that ignore `TREETIME_CARGO_HOME`, because the harness wrapper always reads it
  4. Validate by hand: `dev/smoke` against a baseline commit with an older `dev/docker/run`; confirm in the log that the harness wrapper path runs and that the binary hash matches a build made with the baseline's own wrapper
- O2, baseline tooling: keep the code; replace the "Decision required" block with a statement in the `dev/smoke` usage text that the baseline revision controls its own build container, then delete this issue
