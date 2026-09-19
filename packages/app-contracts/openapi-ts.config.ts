import { defineConfig } from "@hey-api/openapi-ts";

export default defineConfig({
  input: "./openapi.yaml",
  output: "./src/generated",
  plugins: [
    { name: "@hey-api/typescript", requests: false, responses: false },
    { name: "zod", compatibilityVersion: 4, requests: false, responses: false },
  ],
});
