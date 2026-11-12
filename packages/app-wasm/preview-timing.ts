import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { parseArgs } from "node:util";

import { zClockPreview, zDataIndex, zRunRecord, zUploadedInput } from "@neherlab/app-contracts";
import * as z from "zod";

import type createModule from "./pkg/app_wasm.js";

const PACKAGE = resolve("packages/app-wasm/pkg");

const DATA = resolve("data");

const JSON_HEADERS: Header[] = [["content-type", "application/json"]];

const BYTES_HEADERS: Header[] = [["content-type", "application/octet-stream"]];

const PREVIEW_INPUTS = new Set(["tree", "metadata"]);

const moduleSchema = z.object({
  // oxlint-disable-next-line unicorn/no-instanceof-builtins -- the imported module belongs to this realm, so instanceof Function is reliable
  default: z.custom<ModuleFactory>((value) => value instanceof Function),
});

const { repeats } = parseArgs({ options: { repeats: { type: "string", default: "3" } } }).values;

const runs = Number(repeats);

assert.ok(Number.isSafeInteger(runs) && runs > 0, `--repeats ${repeats}`);

const createSingleModule = moduleSchema.parse(await import(join(PACKAGE, "st", "app_wasm.js"))).default;

const wasm = await createSingleModule();

const backend = new wasm.Backend();

const index = zDataIndex.parse(JSON.parse(readFileSync(join(DATA, "index.json"), "utf8")));

const datasets = index.datasets
  .filter((dataset) => dataset.inputs.some((input) => input.kind === "metadata"))
  .toSorted((left, right) => left.samples - right.samples);

writeRow(["dataset", "tips", "median_ms", "min_ms", "max_ms"]);

for (const dataset of datasets) {
  const timing = await timePreview(dataset.inputs.filter((input) => PREVIEW_INPUTS.has(input.kind)));

  writeRow(
    "error" in timing
      ? [dataset.name, String(dataset.samples), timing.error]
      : [dataset.name, String(timing.tips), ...[timing.median, timing.min, timing.max].map((ms) => ms.toFixed(0))],
  );
}

async function timePreview(inputs: { kind: string; path: string }[]): Promise<Timing | Failure> {
  const created = parse(
    zRunRecord,
    await call("POST", "/api/runs", encode({ command: "timetree", config: {}, defer_start: true }), JSON_HEADERS),
  );

  const config: Record<string, string> = {};

  for (const input of inputs) {
    const name = input.path.split("/").at(-1) ?? input.path;

    const response = await call(
      "PUT",
      `/api/runs/${created.id}/inputs/${name}`,
      readFileSync(join(DATA, input.path)),
      BYTES_HEADERS,
    );

    config[input.kind] = parse(zUploadedInput, response).path;
  }

  const body = encode({ command: "timetree", config, beyond_limit: true });
  const durations: number[] = [];
  let tips = 0;

  for (let repeat = 0; repeat < runs; repeat += 1) {
    const start = performance.now();
    const reply = await call("POST", "/api/clock-preview", body, JSON_HEADERS);

    durations.push(performance.now() - start);

    if (reply.status !== 200) {
      return { error: new TextDecoder().decode(reply.body) };
    }

    const preview = parse(zClockPreview, reply);

    tips = preview.tips;
    assert.ok(preview.clock !== undefined, `the preview of ${JSON.stringify(inputs)} has a regression`);
  }

  const sorted = durations.toSorted((left, right) => left - right);

  return {
    tips,
    median: sorted[Math.floor(sorted.length / 2)] ?? Number.NaN,
    min: sorted[0] ?? Number.NaN,
    max: sorted.at(-1) ?? Number.NaN,
  };
}

function writeRow(cells: string[]): void {
  process.stdout.write(`${cells.join("\t")}\n`);
}

async function call(method: string, url: string, body?: Uint8Array, headers: Header[] = []): Promise<Reply> {
  // oxlint-disable-next-line custom/require-io-timeout -- the back end answers in this process, without network IO
  const exchange = backend.fetch(
    method,
    url,
    headers.map(([name]) => name),
    headers.map(([, value]) => value),
    body,
  );

  const head = await exchange.head();
  const chunks: Uint8Array[] = [];

  for (let chunk = await exchange.chunk(); chunk !== undefined; chunk = await exchange.chunk()) {
    chunks.push(chunk);
  }

  exchange.free();

  return { status: head.status, body: Buffer.concat(chunks) };
}

function parse<T extends z.ZodType>(schema: T, reply: Reply): z.infer<T> {
  const text = new TextDecoder().decode(reply.body);

  assert.equal(reply.status, 200, text);

  return schema.parse(JSON.parse(text));
}

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- the request bodies are plain JSON values of the API
function encode(value: unknown): Uint8Array {
  return new TextEncoder().encode(JSON.stringify(value));
}

type ModuleFactory = typeof createModule;

type Header = [name: string, value: string];

interface Reply {
  status: number;
  body: Uint8Array;
}

interface Failure {
  error: string;
}

interface Timing {
  tips: number;
  median: number;
  min: number;
  max: number;
}
