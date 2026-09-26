import { defineConfig } from "@hey-api/openapi-ts";

const LONG_INTEGER_FORMATS = new Set(["int64", "uint64"]);

const BINARY_FORMAT = "binary";

const SAFE_INTEGER_LIMITS = {
  int64: {
    minValue: Number.MIN_SAFE_INTEGER,
    minError: `Invalid value: Expected int64 to be >= ${Number.MIN_SAFE_INTEGER}, the smallest exact JSON integer`,
    maxValue: Number.MAX_SAFE_INTEGER,
    maxError: `Invalid value: Expected int64 to be <= ${Number.MAX_SAFE_INTEGER}, the largest exact JSON integer`,
  },
  uint64: {
    minValue: 0,
    minError: "Invalid value: Expected uint64 to be >= 0",
    maxValue: Number.MAX_SAFE_INTEGER,
    maxError: `Invalid value: Expected uint64 to be <= ${Number.MAX_SAFE_INTEGER}, the largest exact JSON integer`,
  },
};

export default defineConfig({
  input: "./openapi.json",
  output: "./src/generated",
  plugins: [
    { name: "@hey-api/client-fetch" },
    { name: "@hey-api/typescript" },
    { name: "@hey-api/sdk", validator: { request: "zod", response: "zod" } },
    {
      name: "zod",
      compatibilityVersion: 4,
      $resolvers: {
        number: (ctx) => {
          const format = ctx.schema.format;

          if (format === undefined || !LONG_INTEGER_FORMATS.has(format)) {
            return undefined;
          }

          const safeIntegers: typeof ctx.utils = {
            ...ctx.utils,
            shouldCoerceToBigInt: () => false,
            maybeBigInt: (value) => ctx.$.fromValue(value),
            getIntegerLimit: () => (format === "int64" ? SAFE_INTEGER_LIMITS.int64 : SAFE_INTEGER_LIMITS.uint64),
          };

          const numberCtx = { ...ctx, utils: safeIntegers };

          for (const node of [numberCtx.nodes.base, numberCtx.nodes.min, numberCtx.nodes.max]) {
            const chain = node(numberCtx);

            if (chain !== undefined) {
              ctx.chain.current = chain;
            }
          }

          return ctx.chain.current;
        },
        string: (ctx) => {
          if (ctx.schema.format !== BINARY_FORMAT) {
            return undefined;
          }

          ctx.chain.current = ctx.chain.current.attr("instanceof").call(ctx.$("Blob"));

          return ctx.chain.current;
        },
      },
    },
  ],
});
