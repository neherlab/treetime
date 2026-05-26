import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { parse as parseYaml } from "yaml";

const scriptDir = dirname(fileURLToPath(import.meta.url));
const pkgDir = resolve(scriptDir, "..");
const specPath = resolve(pkgDir, "openapi.yaml");
const outDir = resolve(pkgDir, "src", "generated");

interface SchemaProperty {
  type?: string;
  $ref?: string;
  items?: SchemaProperty;
  enum?: string[];
  additionalProperties?: boolean | SchemaProperty;
  properties?: Record<string, SchemaProperty>;
  required?: string[];
}

interface Schema {
  type?: string;
  required?: string[];
  properties?: Record<string, SchemaProperty>;
  items?: SchemaProperty;
  enum?: string[];
  additionalProperties?: boolean | SchemaProperty;
}

interface PathOperation {
  operationId: string;
  "x-bridge-type": "query" | "command";
  requestBody?: {
    content: { "application/json": { schema: SchemaProperty } };
  };
  responses: {
    "200": {
      content: Record<string, { schema: SchemaProperty }>;
    };
  };
}

interface OpenAPISpec {
  paths: Record<string, Record<string, PathOperation>>;
  components: { schemas: Record<string, Schema> };
}

function resolveRef(ref: string): string {
  return ref.replace("#/components/schemas/", "");
}

function schemaToTs(prop: SchemaProperty, required: boolean): string {
  const opt = required ? "" : "?";
  const type = schemaTypeToTs(prop);
  return `${opt}: ${type}`;
}

function schemaTypeToTs(prop: SchemaProperty): string {
  if (prop.$ref) {
    return resolveRef(prop.$ref);
  }
  if (prop.enum) {
    return prop.enum.map((v) => `"${v}"`).join(" | ");
  }
  if (prop.type === "array") {
    const itemType = prop.items ? schemaTypeToTs(prop.items) : "unknown";
    return `${itemType}[]`;
  }
  if (prop.type === "object") {
    if (prop.additionalProperties) {
      return "Record<string, unknown>";
    }
    if (prop.properties) {
      const fields = Object.entries(prop.properties)
        .map(([name, p]) => {
          const req = prop.required?.includes(name) ?? false;
          return `  ${name}${schemaToTs(p, req)};`;
        })
        .join("\n");
      return `{\n${fields}\n}`;
    }
    return "Record<string, never>";
  }
  if (prop.type === "integer" || prop.type === "number") return "number";
  if (prop.type === "boolean") return "boolean";
  if (prop.type === "string") return "string";
  return "unknown";
}

function generateTypes(schemas: Record<string, Schema>): string {
  const lines: string[] = ["// Generated from openapi.yaml - do not edit", ""];

  for (const [name, schema] of Object.entries(schemas)) {
    if (schema.enum) {
      lines.push(`export type ${name} = ${schema.enum.map((v) => `"${v}"`).join(" | ")};`);
      lines.push("");
      continue;
    }

    if (schema.type === "object") {
      lines.push(`export interface ${name} {`);
      if (schema.properties) {
        for (const [propName, prop] of Object.entries(schema.properties)) {
          const req = schema.required?.includes(propName) ?? false;
          lines.push(`  ${propName}${schemaToTs(prop, req)};`);
        }
      }
      lines.push("}");
      lines.push("");
    }
  }

  return lines.join("\n");
}

interface EndpointInfo {
  operationId: string;
  bridgeType: "query" | "command";
  argsType: string | undefined;
  resultType: string;
}

function extractEndpoints(spec: OpenAPISpec): EndpointInfo[] {
  const endpoints: EndpointInfo[] = [];

  for (const methods of Object.values(spec.paths)) {
    for (const op of Object.values(methods)) {
      const argsRef = op.requestBody?.content["application/json"]?.schema.$ref;
      const responseContent = op.responses["200"].content;
      const resultSchema =
        responseContent["application/json"]?.schema ?? responseContent["text/event-stream"]?.schema;
      let resultType: string;
      if (resultSchema?.$ref) {
        resultType = resolveRef(resultSchema.$ref);
      } else if (resultSchema?.type === "array" && resultSchema.items?.$ref) {
        resultType = `${resolveRef(resultSchema.items.$ref)}[]`;
      } else {
        resultType = "unknown";
      }

      endpoints.push({
        operationId: op.operationId,
        bridgeType: op["x-bridge-type"],
        argsType: argsRef ? resolveRef(argsRef) : undefined,
        resultType,
      });
    }
  }

  return endpoints;
}

function generateBridge(endpoints: EndpointInfo[]): string {
  const lines: string[] = ["// Generated from openapi.yaml - do not edit", ""];

  const imports = new Set<string>();
  for (const ep of endpoints) {
    if (ep.argsType) imports.add(ep.argsType);
    if (ep.resultType !== "unknown") imports.add(ep.resultType.replace("[]", ""));
  }
  lines.push(`import type { ${[...imports].sort().join(", ")} } from "./types";`);
  lines.push("");

  lines.push('import type { ProgressEvent } from "./types";');
  lines.push("");

  lines.push("export interface CommandOptions {");
  lines.push("  onProgress?: (event: ProgressEvent) => void;");
  lines.push("  signal?: AbortSignal;");
  lines.push("}");
  lines.push("");

  lines.push("export class CancelledError extends Error {");
  lines.push('  constructor() { super("Operation cancelled"); this.name = "CancelledError"; }');
  lines.push("}");
  lines.push("");

  lines.push("export interface BridgeTransport {");
  lines.push("  query<T>(endpoint: string): Promise<T>;");
  lines.push("  command<T>(endpoint: string, args: unknown, options?: CommandOptions): Promise<T>;");
  lines.push("}");
  lines.push("");

  const queryEndpoints = endpoints.filter((e) => e.bridgeType === "query");
  const commandEndpoints = endpoints.filter((e) => e.bridgeType === "command");

  lines.push("export interface TreeTimeBridge {");
  for (const ep of queryEndpoints) {
    lines.push(`  ${ep.operationId}(): Promise<${ep.resultType}>;`);
  }
  for (const ep of commandEndpoints) {
    lines.push(`  ${ep.operationId}(args: ${ep.argsType}, options?: CommandOptions): Promise<${ep.resultType}>;`);
  }
  lines.push("}");
  lines.push("");

  lines.push("export function createBridge(transport: BridgeTransport): TreeTimeBridge {");
  lines.push("  return {");
  for (const ep of queryEndpoints) {
    lines.push(`    ${ep.operationId}: () => transport.query<${ep.resultType}>("${ep.operationId}"),`);
  }
  for (const ep of commandEndpoints) {
    lines.push(
      `    ${ep.operationId}: (args, options) => transport.command<${ep.resultType}>("${ep.operationId}", args, options),`,
    );
  }
  lines.push("  };");
  lines.push("}");
  lines.push("");

  return lines.join("\n");
}

function generateIndex(schemaNames: string[]): string {
  const lines: string[] = ["// Generated from openapi.yaml - do not edit", ""];
  lines.push('export type { BridgeTransport, CommandOptions, TreeTimeBridge } from "./bridge";');
  lines.push('export { CancelledError, createBridge } from "./bridge";');
  lines.push("");

  lines.push(`export type { ${[...schemaNames].sort().join(", ")} } from "./types";`);
  lines.push("");

  return lines.join("\n");
}

function main() {
  const specText = readFileSync(specPath, "utf-8");
  const spec = parseYaml(specText) as OpenAPISpec;

  mkdirSync(outDir, { recursive: true });

  const endpoints = extractEndpoints(spec);

  writeFileSync(resolve(outDir, "types.ts"), generateTypes(spec.components.schemas));
  writeFileSync(resolve(outDir, "bridge.ts"), generateBridge(endpoints));
  writeFileSync(resolve(outDir, "index.ts"), generateIndex(Object.keys(spec.components.schemas)));

  const queryCount = endpoints.filter((e) => e.bridgeType === "query").length;
  const commandCount = endpoints.filter((e) => e.bridgeType === "command").length;
  console.log(`Generated ${queryCount} queries, ${commandCount} commands, ${Object.keys(spec.components.schemas).length} types`);
}

main();
