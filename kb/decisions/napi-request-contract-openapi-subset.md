# N-API request contract normalized to the openapi subset

The N-API client accepts one request object per operation (`ancestral`, `clock`, `mugration`, `optimize`, `prune`, `timetree`), deserialized from the JSON string passed to the export. Each request struct is the subset of fields the browser client sends, matching the operation schema in `packages/app-contracts/openapi.yaml` and the request struct the HTTP server accepts for the same operation.

## Behavior

Each N-API request struct defines exactly the fields for its operation and maps them to the fixed default output set: the run writes the command's default `--output-all` file set into the request `outdir`, with an empty output selection and the default topology ordering (descendant-count ladderization). Enum-valued fields (alphabet, model, gap-fill policy, ancestral method, branch-length mode, time-marginal mode, reroot method, optimization method, initial-guess mode) deserialize directly into the core domain enums, whose serde spellings are the wire contract. The server and N-API clients accept the same request fields for each operation.

## Alternatives considered

Deserialize the full command-line argument object (the raw CLI parser struct with every flag). Rejected: the raw struct flattens several sub-argument groups, and serde `flatten` disables `deny_unknown_fields`, so a request carrying a field outside the intended UI subset was silently accepted and ignored rather than rejected. The per-operation subset struct rejects unknown fields and documents the exact contract the client depends on.

## Implementation

- `packages/app-napi/src/commands/<op>.rs`: request struct and the `run_<op>` orchestration for each operation
- `packages/app-napi/src/exports.rs`: deserializes the request struct from the argument JSON and runs the orchestration
- `packages/app-server/src/commands/<op>.rs`: the symmetric server request struct and orchestration
