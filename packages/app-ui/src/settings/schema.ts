import { openApiDocument, type AppCommand } from "@neherlab/app-contracts";
import * as z from "zod";

import { isJsonObject, zJsonValue, type JsonValue } from "./json";

export type SettingKind = "switch" | "tristate" | "enum" | "integer" | "number" | "text" | "list" | "enum-list";

type PathRole = "input" | "input-template" | "output";

export type ListItemKind = "string" | "number" | "integer";

export interface SettingOption {
  value: string;
  help: string;
}

export interface SettingSpec {
  key: string;
  path: string[];
  kind: SettingKind;
  nullable: boolean;
  options: SettingOption[];
  itemKind: ListItemKind;
  defaultValue: JsonValue;
  flag: string;
  numArgs: [number, number | null];
  valueDelimiter: string | null;
  cliValues: Record<string, string>;
  pathRole: PathRole | null;
  minimum: number | null;
  help: string;
  more: string;
}

export interface CommandSettings {
  command: AppCommand;
  specs: SettingSpec[];
}

interface SchemaNode {
  $ref?: string | undefined;
  type?: string | string[] | undefined;
  enum?: string[] | undefined;
  const?: string | undefined;
  oneOf?: SchemaNode[] | undefined;
  anyOf?: SchemaNode[] | undefined;
  items?: SchemaNode | undefined;
  properties?: Record<string, SchemaNode> | undefined;
  description?: string | undefined;
  default?: JsonValue | undefined;
  minimum?: number | undefined;
  "x-cli-flag"?: string | undefined;
  "x-cli-num-args"?: [number, number | null] | undefined;
  "x-cli-value-delimiter"?: string | undefined;
  "x-cli-values"?: Record<string, string> | undefined;
  "x-path"?: PathRole | undefined;
}

const zSchemaNode: z.ZodType<SchemaNode> = z.lazy(() =>
  z.object({
    $ref: z.string().optional(),
    type: z.union([z.string(), z.array(z.string())]).optional(),
    enum: z.array(z.string()).optional(),
    const: z.string().optional(),
    oneOf: z.array(zSchemaNode).optional(),
    anyOf: z.array(zSchemaNode).optional(),
    items: zSchemaNode.optional(),
    properties: z.record(z.string(), zSchemaNode).optional(),
    description: z.string().optional(),
    default: zJsonValue.optional(),
    minimum: z.number().optional(),
    "x-cli-flag": z.string().optional(),
    "x-cli-num-args": z.tuple([z.number(), z.number().nullable()]).optional(),
    "x-cli-value-delimiter": z.string().optional(),
    "x-cli-values": z.record(z.string(), z.string()).optional(),
    "x-path": z.enum(["input", "input-template", "output"]).optional(),
  }),
);

const zComponents = z.object({ components: z.object({ schemas: z.record(z.string(), zJsonValue) }) });

const SCHEMA_KEY = "$schema";

const CONFIG_COMPONENTS: Record<AppCommand, string> = {
  timetree: "TimetreeConfig",
  optimize: "OptimizeConfig",
  prune: "PruneConfig",
  ancestral: "AncestralConfig",
  clock: "ClockConfig",
  mugration: "MugrationConfig",
};

class SchemaError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "SchemaError";
  }
}

export function commandSettings(command: AppCommand): CommandSettings {
  const components = zComponents.parse(openApiDocument).components.schemas;

  const resolve = (name: string): SchemaNode => {
    const node = components[name];

    if (node === undefined) {
      throw new SchemaError(`the OpenAPI document has no component \`${name}\``);
    }

    return zSchemaNode.parse(node);
  };

  const root = resolve(CONFIG_COMPONENTS[command]);

  return { command, specs: collectSpecs(root, [], null, resolve) };
}

function collectSpecs(
  node: SchemaNode,
  parentPath: string[],
  parentDefault: JsonValue,
  resolve: (name: string) => SchemaNode,
): SettingSpec[] {
  return Object.entries(node.properties ?? {}).flatMap(([key, property]) => {
    if (parentPath.length === 0 && key === SCHEMA_KEY) {
      return [];
    }

    const path = [...parentPath, key];
    const inheritedDefault = childDefault(parentDefault, key);
    const nested = nestedObject(property, resolve);

    if (nested !== undefined) {
      return collectSpecs(nested, path, property.default ?? inheritedDefault, resolve);
    }

    return [settingSpec(property, path, inheritedDefault, resolve)];
  });
}

function childDefault(parentDefault: JsonValue, key: string): JsonValue {
  return isJsonObject(parentDefault) ? (parentDefault[key] ?? null) : null;
}

function nestedObject(property: SchemaNode, resolve: (name: string) => SchemaNode): SchemaNode | undefined {
  const candidates = [property, ...(property.anyOf ?? []), ...(property.oneOf ?? [])];

  for (const candidate of candidates) {
    if (candidate.$ref !== undefined) {
      const target = resolve(refName(candidate.$ref));

      if (target.properties !== undefined) {
        return target;
      }
    }
  }

  return undefined;
}

function settingSpec(
  property: SchemaNode,
  path: string[],
  inheritedDefault: JsonValue,
  resolve: (name: string) => SchemaNode,
): SettingSpec {
  const key = path.join(".");
  const flag = property["x-cli-flag"];
  const numArgs = property["x-cli-num-args"];

  if (flag === undefined || numArgs === undefined) {
    throw new SchemaError(`setting \`${key}\` has no command-line annotation`);
  }

  const { help, more } = splitDescription(property.description ?? "");
  const nullable = isNullable(property);
  const base = nonNullBranch(property, resolve);
  const options = enumOptions(base, resolve);
  const types = typeList(base);

  const spec: SettingSpec = {
    key,
    path,
    kind: "text",
    nullable,
    options: [],
    itemKind: "string",
    defaultValue: property.default ?? inheritedDefault,
    flag,
    numArgs,
    valueDelimiter: property["x-cli-value-delimiter"] ?? null,
    cliValues: property["x-cli-values"] ?? {},
    pathRole: property["x-path"] ?? null,
    minimum: base.minimum ?? null,
    help,
    more,
  };

  if (options.length > 0) {
    return { ...spec, kind: "enum", options };
  }

  if (types.includes("array")) {
    const items = base.items === undefined ? {} : deref(base.items, resolve);
    const itemOptions = enumOptions(items, resolve);

    if (itemOptions.length > 0) {
      return { ...spec, kind: "enum-list", options: itemOptions };
    }

    return { ...spec, kind: "list", itemKind: listItemKind(typeList(items), key) };
  }

  if (types.includes("boolean")) {
    return { ...spec, kind: nullable ? "tristate" : "switch" };
  }

  if (types.includes("integer")) {
    return { ...spec, kind: "integer" };
  }

  if (types.includes("number")) {
    return { ...spec, kind: "number" };
  }

  if (types.includes("string")) {
    return spec;
  }

  throw new SchemaError(`setting \`${key}\` has a schema construct the form does not render`);
}

function listItemKind(types: string[], key: string): ListItemKind {
  if (types.includes("integer")) {
    return "integer";
  }

  if (types.includes("number")) {
    return "number";
  }

  if (types.includes("string")) {
    return "string";
  }

  throw new SchemaError(`list setting \`${key}\` has items the form does not render`);
}

function isNullable(property: SchemaNode): boolean {
  return (
    typeList(property).includes("null") ||
    [...(property.anyOf ?? []), ...(property.oneOf ?? [])].some((branch) => typeList(branch).includes("null"))
  );
}

function nonNullBranch(property: SchemaNode, resolve: (name: string) => SchemaNode): SchemaNode {
  const branches = property.anyOf?.filter((branch) => !typeList(branch).includes("null"));

  if (branches !== undefined && branches.length === 1 && branches[0] !== undefined) {
    return deref(branches[0], resolve);
  }

  return deref(property, resolve);
}

function enumOptions(node: SchemaNode, resolve: (name: string) => SchemaNode): SettingOption[] {
  const target = deref(node, resolve);
  const direct = (target.enum ?? []).map((value) => ({ value, help: "" }));

  const constant =
    target.const === undefined ? [] : [{ value: target.const, help: splitDescription(target.description ?? "").help }];

  const alternatives = (target.oneOf ?? []).flatMap((branch) => enumOptions(branch, resolve));

  return [...direct, ...constant, ...alternatives];
}

function deref(node: SchemaNode, resolve: (name: string) => SchemaNode): SchemaNode {
  return node.$ref === undefined ? node : resolve(refName(node.$ref));
}

function typeList(node: SchemaNode): string[] {
  if (node.type === undefined) {
    return [];
  }

  return Array.isArray(node.type) ? node.type : [node.type];
}

function refName(reference: string): string {
  return reference.slice(reference.lastIndexOf("/") + 1);
}

interface DescriptionParts {
  help: string;
  more: string;
}

function splitDescription(description: string): DescriptionParts {
  const paragraphs = description.split(/\n\s*\n/u).flatMap((paragraph) => {
    const joined = paragraph.replaceAll(/\s*\n\s*/gu, " ").trim();

    return joined === "" ? [] : [joined];
  });

  return { help: paragraphs[0] ?? "", more: paragraphs.slice(1).join("\n\n") };
}
