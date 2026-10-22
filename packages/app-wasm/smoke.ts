import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { parseArgs } from "node:util";

import { zRunList, zRunRecord, zRunsFilesResponse, zUploadedInput, zVersionInfo } from "@neherlab/app-contracts";
import * as z from "zod";

import type createModule from "./pkg/app_wasm.js";

const PACKAGE = resolve("packages/app-wasm/pkg");

const DATASET = resolve("data/flu/h3n2/20");

const JSON_HEADERS: Header[] = [["content-type", "application/json"]];

const BYTES_HEADERS: Header[] = [["content-type", "application/octet-stream"]];

const ZIP_SIGNATURE = "PK";

const KEEP_ALIVE_INTERVAL_MS = 1000;

const moduleSchema = z.object({
  // oxlint-disable-next-line unicorn/no-instanceof-builtins -- the imported module belongs to this realm, so instanceof Function is reliable
  default: z.custom<ModuleFactory>((value) => value instanceof Function),
});

const { jobs } = parseArgs({ options: { jobs: { type: "string" } } }).values;

const threads = jobs === undefined ? undefined : Number(jobs);

assert.ok(threads === undefined || (Number.isSafeInteger(threads) && threads > 0), `--jobs ${String(jobs)}`);

const createSingleModule = await moduleFactory("st");

const wasm = await createSingleModule();

const singleJobWasm = threads === undefined ? await createSingleModule() : undefined;

const backend = new wasm.Backend();

await checkNewBackend();

await checkAncestralRun();

await checkCancelledRun();

await checkChangedInput();

await checkEmptyInput();

await checkRepeatedInput();

await checkClockRunWithoutCharts();

await checkEventStreamAborts();

async function checkNewBackend(): Promise<void> {
  const version = await call("GET", "/api/version");
  const runs = await call("GET", "/api/runs");

  assert.deepEqual(
    [version.status, runs.status, parse(zRunList, runs).runs],
    [200, 200, []],
    `${text(version)} ${text(runs)}`,
  );
  parse(zVersionInfo, version);
}

async function checkAncestralRun(): Promise<void> {
  const created = parse(
    zRunRecord,
    await call("POST", "/api/runs", encode({ command: "ancestral", config: {}, defer_start: true }), JSON_HEADERS),
  );

  const tree = await upload(created.id, "tree.nwk");
  const alignment = await upload(created.id, "aln.fasta.xz");

  const started = await call(
    "POST",
    `/api/runs/${created.id}/start`,
    encode({ config: { tree, alignment: [alignment] } }),
    JSON_HEADERS,
  );

  const computing = parse(zRunRecord, await call("GET", `/api/runs/${created.id}`));
  const [order] = backend.takeJobOrders();
  assert.ok(order?.job !== undefined, "the start of a run queues a start order");
  const inputPaths = order.inputPaths;
  const outputs = await driveJob(created.id, order.job, inputPaths);

  const record = parse(zRunRecord, await call("GET", `/api/runs/${created.id}`));
  const files = parse(zRunsFilesResponse, await call("GET", `/api/runs/${created.id}/files`));
  const newick = await call("GET", `/api/runs/${created.id}/file?path=ancestral.nwk`);
  const archive = await call("GET", `/api/runs/${created.id}/archive`);

  assert.deepEqual(
    [
      started.status,
      computing.status,
      inputPaths,
      outputs.completed,
      record.status,
      files.map((file) => file.path),
      record.output_files.every((file) => outputs.files.includes(file.path)),
      files.some((file) => file.path === "ancestral.sequences.fasta"),
      newick.status,
      text(newick).slice(0, 1),
      archive.status,
      text(archive).slice(0, ZIP_SIGNATURE.length),
    ],
    [200, "running", [alignment, tree], true, "ok", outputs.files.toSorted(), true, true, 200, "(", 200, ZIP_SIGNATURE],
    text(started),
  );
}

async function checkCancelledRun(): Promise<void> {
  const id = await startedRun("ancestral", async (runId) => ({
    tree: await upload(runId, "tree.nwk"),
    alignment: [await upload(runId, "aln.fasta.xz")],
  }));

  const [start] = backend.takeJobOrders();
  const cancel = await call("POST", `/api/runs/${id}/cancel`);
  const orders = backend.takeJobOrders().map((order) => [order.kind, order.runId]);
  const failed = backend.failRun(id, []);

  const completed = backend.completeRun(
    id,
    JSON.stringify({ terminal: { status: "cancelled", job_id: id }, warnings: [] }),
  );

  const record = parse(zRunRecord, await call("GET", `/api/runs/${id}`));

  assert.deepEqual(
    [start?.kind, cancel.status, orders, failed, completed, record.status],
    ["start", 200, [["cancel", id]], true, false, "cancelled"],
  );
}

async function checkChangedInput(): Promise<void> {
  const id = await startedRun("ancestral", async (runId) => ({
    tree: await upload(runId, "tree.nwk"),
    alignment: [await upload(runId, "aln.fasta.xz")],
  }));

  const [order] = backend.takeJobOrders();
  assert.ok(order?.job !== undefined, "the start of a run queues a start order");
  await driveJob(id, order.job, order.inputPaths, (path, bytes) =>
    path.endsWith("tree.nwk") ? Buffer.concat([bytes, encode(" ")]) : bytes,
  );
  const record = parse(zRunRecord, await call("GET", `/api/runs/${id}`));

  assert.equal(record.status, "error");
  assert.match(record.error?.message ?? "", /^input '.+tree\.nwk' differs from the file the run prepared/u);
}

async function checkEmptyInput(): Promise<void> {
  const id = await startedRun("prune", async (runId) => ({
    tree: await upload(runId, "tree.nwk"),
    prune_nodes_list_file: await uploadBytes(runId, "empty.txt", new Uint8Array()),
  }));

  const [order] = backend.takeJobOrders();
  assert.ok(order?.job !== undefined, "the start of a run queues a start order");
  const outputs = await driveJob(id, order.job, order.inputPaths);
  const record = parse(zRunRecord, await call("GET", `/api/runs/${id}`));

  assert.deepEqual([outputs.completed, record.status], [true, "ok"], JSON.stringify(record.error));
}

async function checkRepeatedInput(): Promise<void> {
  const uploads: string[] = [];

  const id = await startedRun("ancestral", async (runId) => {
    const alignment = await upload(runId, "aln.fasta.xz");
    const tree = await upload(runId, "tree.nwk");

    uploads.push(alignment, tree);

    return { tree, alignment: [alignment, alignment] };
  });

  const [order] = backend.takeJobOrders();
  assert.ok(order?.job !== undefined, "the start of a run queues a start order");
  const inputPaths = order.inputPaths;
  const outputs = await driveJob(id, order.job, inputPaths);
  const record = parse(zRunRecord, await call("GET", `/api/runs/${id}`));

  assert.deepEqual([inputPaths, outputs.completed, record.status], [uploads, true, "ok"], JSON.stringify(record.error));
}

async function checkClockRunWithoutCharts(): Promise<void> {
  const id = await startedRun("clock", async (runId) => ({
    tree: await upload(runId, "tree.nwk"),
    metadata: await upload(runId, "metadata.tsv"),
  }));

  const [order] = backend.takeJobOrders();
  assert.ok(order?.job !== undefined, "the start of a clock run queues a start order");
  const outputs = await driveJob(id, order.job, order.inputPaths);

  const record = parse(zRunRecord, await call("GET", `/api/runs/${id}`));
  const files = parse(zRunsFilesResponse, await call("GET", `/api/runs/${id}/files`));
  const paths = new Set(files.map((file) => file.path));

  assert.deepEqual(
    [
      outputs.completed,
      record.status,
      ["clock.clock-model.json", "clock.clock.csv", "clock.nwk"].filter((path) => !paths.has(path)),
      ["clock.svg", "clock.png", "clock.pdf"].filter((path) => paths.has(path)),
    ],
    [true, "ok", [], []],
    JSON.stringify(record.error),
  );
}

async function startedRun(command: string, config: (id: string) => Promise<RunInputs>): Promise<string> {
  const created = parse(
    zRunRecord,
    await call("POST", "/api/runs", encode({ command, config: {}, defer_start: true }), JSON_HEADERS),
  );

  const started = await call(
    "POST",
    `/api/runs/${created.id}/start`,
    encode({ config: await config(created.id) }),
    JSON_HEADERS,
  );

  assert.equal(started.status, 200, text(started));

  return created.id;
}

async function driveJob(
  id: string,
  jobJson: string,
  inputPaths: string[],
  change: (path: string, bytes: Uint8Array) => Uint8Array = (_path, bytes) => bytes,
): Promise<{ completed: boolean; files: string[] }> {
  const jobWasm = singleJobWasm ?? (await threadedJobModule());
  const job = new jobWasm.Job(jobJson, threads ?? 1);

  for (const path of inputPaths) {
    const reader = backend.readRunInput(id, path);
    const bytes = change(path, Buffer.concat(readAll(reader)));

    reader.free();
    job.beginInput(path);

    if (bytes.length > 0) {
      job.appendInput(path, bytes);
    }
  }

  const outcome = await keepingNodeAlive(
    job.run((event) => {
      backend.appendRunEvent(id, event);
    }),
  );

  const files = job.outputFiles();

  for (const path of files) {
    const reader = job.readOutput(path);

    backend.beginRunOutput(id, path);
    readAll(reader).forEach((chunk) => {
      backend.appendRunOutput(id, path, chunk);
    });
    reader.free();
  }

  job.free();

  return { completed: backend.completeRun(id, outcome), files };
}

async function keepingNodeAlive<T>(pending: Promise<T>): Promise<T> {
  const timer = setInterval(() => undefined, KEEP_ALIVE_INTERVAL_MS);

  try {
    return await pending;
  } finally {
    clearInterval(timer);
  }
}

async function threadedJobModule(): Promise<JobModule> {
  const createThreadedModule = await moduleFactory("mt");

  Object.assign(globalThis, { treetimeJobs: threads });

  return createThreadedModule();
}

async function moduleFactory(build: "st" | "mt"): Promise<ModuleFactory> {
  return moduleSchema.parse(await import(join(PACKAGE, build, "app_wasm.js"))).default;
}

function readAll(reader: { chunk(): Uint8Array | undefined }): Uint8Array[] {
  const chunks: Uint8Array[] = [];

  for (let chunk = reader.chunk(); chunk !== undefined; chunk = reader.chunk()) {
    chunks.push(chunk);
  }

  return chunks;
}

async function checkEventStreamAborts(): Promise<void> {
  // oxlint-disable-next-line custom/require-io-timeout -- the back end answers in this process; the check aborts the stream itself
  const exchange = backend.fetch("GET", "/api/events", [], [], undefined);
  const head = await exchange.head();
  const contentType = head.headerValues[head.headerNames.indexOf("content-type")];

  exchange.abort();
  exchange.free();
  assert.deepEqual([head.status, contentType], [200, "text/event-stream"]);
}

async function upload(id: string, name: string): Promise<string> {
  return uploadBytes(id, name, readFileSync(join(DATASET, name)));
}

async function uploadBytes(id: string, name: string, bytes: Uint8Array): Promise<string> {
  const response = await call("PUT", `/api/runs/${id}/inputs/${name}`, bytes, BYTES_HEADERS);

  return parse(zUploadedInput, response).path;
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
  return schema.parse(JSON.parse(text(reply)));
}

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- the request bodies are plain JSON values of the API
function encode(value: unknown): Uint8Array {
  return new TextEncoder().encode(JSON.stringify(value));
}

function text(reply: Reply): string {
  return new TextDecoder().decode(reply.body);
}

type ModuleFactory = typeof createModule;

type JobModule = Awaited<ReturnType<ModuleFactory>>;

type Header = [name: string, value: string];

type RunInputs = Record<string, string | string[]>;

interface Reply {
  status: number;
  body: Uint8Array;
}
