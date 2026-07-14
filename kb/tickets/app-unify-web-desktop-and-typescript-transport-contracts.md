# Unify web, desktop, and TypeScript transport contracts

Use one generated request/result contract and one Rust conversion layer for every application client.

## Required changes

- Add all optimize reroot fields and serialized value constraints to OpenAPI.
- Generate TypeScript DTOs from the updated schema.
- Define a discriminated success/error/cancelled envelope; reject failed promises in the web bridge.
- Deserialize the same transport DTO in desktop and web, then call one explicit DTO-to-core conversion.
- Reject unknown fields at transport boundaries where they indicate schema drift.

## Validation

- Schema-generation consistency check.
- Web and desktop parity for every command and argument field.
- Success, validation error, panic, cancellation, and unknown-field cases.
- Rust, TypeScript, lint, formatting, and test suites.

## Related issues

- Source: [kb/issues/H-app-transport-contracts-diverge-across-clients.md](../issues/H-app-transport-contracts-diverge-across-clients.md)
