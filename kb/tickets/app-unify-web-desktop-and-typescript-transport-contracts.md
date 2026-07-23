# Unify web, desktop, and TypeScript transport contracts

Use one generated transport contract and one Rust conversion layer for every application client.

## Required changes

- Define canonical request, result, progress-event, error, and cancellation DTOs for every command.
- Add all optimize reroot fields and serialized value constraints to the schema.
- Generate TypeScript DTOs from the canonical schema and fail validation when generated output is stale.
- Use one discriminated success/error/cancelled envelope across HTTP streaming and NAPI promises.
- Deserialize the same transport DTO in desktop and web, then call one explicit DTO-to-application conversion.
- Reject unknown fields at transport boundaries where they indicate schema drift.
- Replace process-global NAPI cancellation with task-scoped cancellation represented by the canonical contract.

## Validation

- Schema-generation consistency check.
- Web and desktop parity for every command and argument field.
- Success, validation error, panic, cancellation, and unknown-field cases.
- Rust, TypeScript, lint, formatting, and test suites.

## Related issues

- Source: [kb/issues/H-app-transport-contracts-diverge-across-clients.md](../issues/H-app-transport-contracts-diverge-across-clients.md)
- Source: [kb/issues/H-app-napi-cancellation-is-process-global.md](../issues/H-app-napi-cancellation-is-process-global.md)
