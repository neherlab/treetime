# Application transport contracts diverge across clients

Web, desktop, Rust, OpenAPI, and generated TypeScript do not share one request/result contract.

## Evidence

- `struct ServerOptimizeArgs` contains reroot fields that are absent from the OpenAPI `OptimizeArgs` schema [packages/app-server/src/args.rs#L397-L456](../../packages/app-server/src/args.rs#L397-L456) [packages/app-contracts/openapi.yaml#L380-L409](../../packages/app-contracts/openapi.yaml#L380-L409).
- `fn command_sse()` emits command failures as `result` events [packages/app-server/src/sse.rs#L96-L120](../../packages/app-server/src/sse.rs#L96-L120), while `function createWebBridge()` resolves every result event [packages/app-web/src/bridge-web.ts#L40-L57](../../packages/app-web/src/bridge-web.ts#L40-L57).
- Desktop forwards JSON strings directly to command-specific N-API deserialization [packages/app-desktop/src/preload.ts#L29-L35](../../packages/app-desktop/src/preload.ts#L29-L35) [packages/app-napi/src/commands.rs#L76-L85](../../packages/app-napi/src/commands.rs#L76-L85).

## Failures

- Optimize reroot fields exist in `ServerOptimizeArgs` but not in OpenAPI or generated `OptimizeArgs`.
- Web SSE sends command failures as `result` events; the TypeScript bridge resolves them as successful typed values.
- Desktop serializes flat transport DTOs directly into nested core command structs. `#[serde(default)]` can hide unknown/dropped fields.

## Potential solutions

- O1. One transport DTO/schema and one explicit conversion shared by web and desktop.
- O2. Separate client DTOs with independent converters and cross-client contract tests. This permits intentional differences but multiplies schema ownership.

## Recommendation

Define one canonical transport DTO and discriminated result envelope. Generate TypeScript from the same schema, deserialize the same DTO for web and desktop, and share one explicit conversion into core command arguments. Core structs are not wire formats.

## Related issues

- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
