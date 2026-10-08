import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { join } from "node:path";

import { zRunList, zRunRecord, zRunsFilesResponse, zUploadedInput, zVersionInfo } from "@neherlab/app-contracts";
import type * as z from "zod";

import createBackendModule from "./pkg/app_wasm.js";

const DATASET = join(import.meta.dirname, "../../data/flu/h3n2/20");

const JSON_HEADERS: Header[] = [["content-type", "application/json"]];

const BYTES_HEADERS: Header[] = [["content-type", "application/octet-stream"]];

const ZIP_SIGNATURE = "PK";

const wasm = await createBackendModule();

wasm.setPanicSink((text) => {
  throw new Error(`the back end panicked: ${text}`);
});

const backend = new wasm.Backend();

await checkNewBackend();

await checkAncestralRun();

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

  const record = parse(zRunRecord, await call("GET", `/api/runs/${created.id}`));
  const files = parse(zRunsFilesResponse, await call("GET", `/api/runs/${created.id}/files`));
  const newick = await call("GET", `/api/runs/${created.id}/file?path=ancestral.nwk`);
  const archive = await call("GET", `/api/runs/${created.id}/archive`);

  assert.deepEqual(
    [
      started.status,
      record.status,
      files.some((file) => file.path === "ancestral.reconstructed-nuc.fasta"),
      newick.status,
      text(newick).slice(0, 1),
      archive.status,
      text(archive).slice(0, ZIP_SIGNATURE.length),
    ],
    [200, "ok", true, 200, "(", 200, ZIP_SIGNATURE],
    text(started),
  );
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
  const response = await call(
    "PUT",
    `/api/runs/${id}/inputs/${name}`,
    readFileSync(join(DATASET, name)),
    BYTES_HEADERS,
  );

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

type Header = [name: string, value: string];

interface Reply {
  status: number;
  body: Uint8Array;
}
