# A malformed `Last-Event-ID` header is silently ignored

Both event streams, `GET /api/runs/{id}/events` and `GET /api/events`, resume after the event named by the `Last-Event-ID` request header. When the header is not a number, the server ignores it without an error and falls back to the `from` query parameter, or to its default when `from` is absent. A malformed `from` parameter, in contrast, is rejected with `400` and an `invalid_request` `ErrorResponse`.

## Evidence

- `fn resume_from()` [packages/app-server/src/routes.rs#L512-L519](../../packages/app-server/src/routes.rs#L512-L519) reads the header with `value.to_str().ok()` and `value.parse::<u64>().ok()`, so a header that is not ASCII or not a number becomes `None`, and then takes `query.from`. The `result_defaulted` lint suppression on the function states this fallback as intended
- `fn runs_events()` uses `resume_from(...).unwrap_or(0)`: the run stream restarts at its first event
- `fn events()` passes `None` on: the app-wide stream sends only events from now on, without a `resync` event

## Impact

- A client that sends a malformed header gets no error. On a run stream it receives every event again; on the app-wide stream it misses the changes between its last event and the reconnect, and its cached answers stay stale
- The app's own client is not affected: `resumableStream` in `packages/app-contracts/src/client.ts` resumes with the `from` query parameter, and the generated SSE client sends `Last-Event-ID` only with the `id` of an event it received

## Possible change

> [!IMPORTANT]
> **Decision required.** The code treats a malformed header as an absent one on purpose (the lint suppression reason on `resume_from()`), while the `from` parameter is rejected. Options:
>
> - Reject a `Last-Event-ID` header that is not a sequence number with `400` and an `invalid_request` `ErrorResponse`, as the `from` parameter is
> - Keep the fallback and treat the header as a hint, as an absent header is
>
> Either way, document the header on both operations and add route tests for a malformed header on both streams.

## Related

- [kb/decisions/app-request-contract-full-command-config.md](../decisions/app-request-contract-full-command-config.md): event streams, stale paths and the UI data layer

## Proposed implementation plan

The cited lines of `fn resume_from()` and its callers (`fn runs_events()` at [packages/app-server/src/routes.rs#L494](../../packages/app-server/src/routes.rs#L494), `fn events()` at [packages/app-server/src/routes.rs#L502](../../packages/app-server/src/routes.rs#L502)) are current. The plan stops at the open decision above; one short plan per option:

- **Reject**: make `fn resume_from()` return a `Result` with the error that a malformed `from` query parameter already produces (find its producer from the `invalid_request` cases of `packages/app-server/src/__tests__/test_routes.rs`), map a header that is not ASCII or not a `u64` to `400` `invalid_request`, propagate with `?` in both handlers, and delete the `result_defaulted` expectation, because nothing is defaulted then
- **Keep the fallback**: no code change; keep the expectation and its reason

For both options: document `Last-Event-ID` on both operations (the route doc comments that feed `packages/app-contracts/openapi.json`, then `just gen`), and add `rstest` route tests in `packages/app-server/src/__tests__/` for a malformed header (`abc`) and a valid one on both streams, asserting the whole response (status and `ErrorResponse`) or the first event `seq`. Oracle: the chosen contract. See red with a wrong expected status.

### Plan validation

- `just lint-rs`, `just test-rs -E 'package(app-server)'`, and `just generated-check` pass
