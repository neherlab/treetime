# The concurrent-upload test does not force the order that the run lock prevents

`fn RunManager.upload_input()` holds the input lock of the run while it computes the remaining upload budget and writes the file [packages/app-backend/src/runs/manager.rs#L223-L240](../../packages/app-backend/src/runs/manager.rs#L223-L240), so a second upload to the same run waits for the first and counts its size. Without the lock, two uploads could both compute their budget before either file is written, and together exceed the input limit of the run.

`test_manager_second_upload_waits_for_the_first_and_counts_it_toward_the_limit` [packages/app-backend/src/runs/__tests__/test_manager.rs#L396-L422](../../packages/app-backend/src/runs/__tests__/test_manager.rs#L396-L422) blocks the first upload inside its body read, starts the second upload in another thread, and releases the first right after it spawns the second. The test cannot tell whether the second thread reached its budget check before the first upload resumed, so on a scheduling where the second thread starts late, the test also passes without the lock. It does not prove that the lock prevents the race.

## Fix direction

> [!IMPORTANT]
> **Investigation required.** A deterministic way to make the second upload reach its budget check while the first is blocked, without code in production that exists only for tests. The test cannot observe a thread that waits on the lock. Candidates: a model-checking test of the lock and budget logic, for example with `loom`, which explores every interleaving; or a budget computation moved into a function that a test can call with both uploads in flight.

## Validation

- With the lock removed locally, the new test fails on every run
- With the lock in place, it passes

## Proposed implementation plan

`fn RunManager.upload_input()` takes `inputs.lock` at [packages/app-backend/src/runs/manager.rs#L235](../../packages/app-backend/src/runs/manager.rs#L235) and holds it through `fn store_upload()`; the test still releases the first upload right after it spawns the second ([packages/app-backend/src/runs/__tests__/test_manager.rs#L405-L406](../../packages/app-backend/src/runs/__tests__/test_manager.rs#L405-L406)). The plan stops at the open investigation above; one short plan per candidate:

- **`loom` model check**: add `loom` through "Adding a dependency" as a dev-dependency of `app-backend`; extract the budget check and the write of `fn store_upload()` behind a small store trait, so a loom test runs two uploads against an in-memory store under `loom::model` and asserts that the total never exceeds the limit. Oracle: the limit contract (sum of stored sizes ≤ `max_total_size`). See red: remove the `inputs.lock.lock()` line and confirm that loom reports the interleaving
- **Testable budget function**: move the budget computation into a pure function of (stored sizes, in-flight size, limit), and test it with both uploads in flight by an `rstest` table. This tests the arithmetic only; the lock itself then stays covered only by the current timing-dependent test, so this option does not meet the validation below on its own

The candidate needs the user's choice after the investigation.
