# The generated zod schemas type optional fields as `T | undefined`

The zod schemas that `@hey-api/openapi-ts` generates into `packages/app-contracts/src/generated/zod.gen.ts` declare each optional field with `.optional()`, whose output type is `T | undefined`. The project compiles with `exactOptionalPropertyTypes: true` (`tsconfig.base.json`), under which the contract types of `types.gen.ts` declare an optional field as absent or `T`, never `undefined`. The parse result therefore does not fit the contract type (the generated folder itself compiles with the option off, `packages/app-contracts/tsconfig.generated.json`).

The consumers work around it: `assertDataIndex()` validates the data index with the schema but returns the unparsed value as the contract type, in [packages/app-ui/build/example-data-vite.ts#L46](../../packages/app-ui/build/example-data-vite.ts#L46) and [packages/app-ui/src/analysis/exampleData.ts#L40](../../packages/app-ui/src/analysis/exampleData.ts#L40). The returned index keeps unknown keys, and no default of the schema is applied, so the value that the app uses can differ from the value that the schema describes.

## Fix direction

Generate zod schemas whose optional fields fit exact optional property types (a generator option or a post-processing step of `just gen` in the zod plugin configuration of `packages/app-contracts/openapi-ts.config.ts`), then return the parse result from `assertDataIndex()` and the other validators instead of the input value.

## Validation

- `zDataIndex.parse()` returns a value of type `DataIndex` under `exactOptionalPropertyTypes`, and both callers use the parsed value
