# Optional schema fields are never nullable

An `Option` field of a type that derives `JsonSchema` is optional and never nullable. This holds in the OpenAPI document, in the command config schemas, and in the JSON schemas of `packages/schemas`.

## Behavior

- **Schemas**: a project schema transform removes `null` from every generated schema: from `type` arrays, from `anyOf` branches, and from `enum` values, and it removes `default: null`. It marks each setting it changed with `x-unset: true`, so the UI knows that the setting can be unset
- **Serialization**: every type that derives both `Serialize` and `JsonSchema` and has `Option` fields leaves out `None` fields, through `#[skip_serializing_none]` of `serde_with` on the type. serde has no global setting for this
- **Input**: a config file or a request that writes `key: null` is rejected by the schema check, with the location of the key
- **UI**: the UI unsets a setting by deleting its key. The generated TypeScript types have `?:` without `null`, and the generated zod schemas use `.optional()`. Under `exactOptionalPropertyTypes`, the output type of `.optional()` still admits `undefined`, so a module that validates untyped data at a trust boundary (stored preferences, the bundled setting catalog) checks it with the generated schema in an assertion function and keeps the generated type, instead of using the zod output type
- **Guards**: tests over the OpenAPI document fail on any schema that allows `null`, on `default: null`, and on an `Option` field that the serialize contract still writes as `null`

## Reason

One rule for missing values: no code converts between `null` and `undefined`, and the generated types state exactly which fields can be missing.

## Implementation

- `packages/treetime-schema/src/no_null.rs`: the schema transform
- `packages/app-server/src/__tests__/test_openapi_contract.rs`: the guards over the OpenAPI document, including a second document built with the serialize contract
- `packages/app-cli/src/cli/schema.rs`: the `default: null` guard over `packages/schemas`
