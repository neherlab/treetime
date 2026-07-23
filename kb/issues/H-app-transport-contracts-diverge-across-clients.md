# Application transport contracts diverge across clients

Web, desktop, Rust, OpenAPI, and generated TypeScript do not share one request, event, and result contract.

## Evidence

- `struct ServerOptimizeArgs` contains reroot fields that are absent from the OpenAPI `OptimizeArgs` schema [`packages/app-server/src/args.rs#L397`](../../packages/app-server/src/args.rs#L397) [`packages/app-contracts/openapi.yaml#L380`](../../packages/app-contracts/openapi.yaml#L380).
- `fn command_sse()` emits command failures as `result` events, while `function createWebBridge()` resolves every result event as the promised result type [`packages/app-server/src/sse.rs#L93`](../../packages/app-server/src/sse.rs#L93) [`packages/app-web/src/bridge-web.ts#L40`](../../packages/app-web/src/bridge-web.ts#L40).
- Desktop forwards JSON strings directly to command-specific N-API deserialization [`packages/app-desktop/src/preload.ts#L29`](../../packages/app-desktop/src/preload.ts#L29) [`packages/app-napi/src/commands.rs#L76`](../../packages/app-napi/src/commands.rs#L76).
- Six server DTOs manually reconstruct nested command argument structs, duplicating fields and defaults from core and OpenAPI [`packages/app-server/src/args.rs#L21`](../../packages/app-server/src/args.rs#L21).
- Schema generation failure is reduced to a Cargo warning, allowing a successful build to retain stale contracts [`packages/app-cli/build.rs#L11`](../../packages/app-cli/build.rs#L11).
- Command defaults and option sets are repeated in core arguments, server arguments, and UI metadata, including timetree `max_iter` [`packages/treetime/src/commands/timetree/args.rs#L111`](../../packages/treetime/src/commands/timetree/args.rs#L111) [`packages/app-server/src/args.rs#L203`](../../packages/app-server/src/args.rs#L203) [`packages/app-ui/src/components/ParamForm.tsx#L15`](../../packages/app-ui/src/components/ParamForm.tsx#L15).

## Failures

- Optimize reroot fields exist in the server DTO but not in OpenAPI or generated TypeScript.
- Web SSE sends command failures through the success event name.
- Core serde defaults can accept absent fields while transport owners silently drift.
- A schema-generation error does not fail the build that consumes generated contracts.

## Required contract

One canonical transport model must own command names, request DTOs, result DTOs, progress/log events, cancellation, and a discriminated success/error/cancelled terminal envelope. HTTP/SSE and Electron/N-API remain transport adapters. Conversion from transport requests to application requests occurs once and rejects unknown or unrepresented fields.

## Related issues

- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
- [H-app-napi-cancellation-is-process-global.md](H-app-napi-cancellation-is-process-global.md)

## Related tickets

- [kb/tickets/app-unify-web-desktop-and-typescript-transport-contracts.md](../tickets/app-unify-web-desktop-and-typescript-transport-contracts.md)
