import { GlobalRegistrator } from "@happy-dom/global-registrator";
import type {
  AppCommand,
  CheckConfigResponse,
  DataIndex,
  RunRecord,
  SparseConfig,
  UploadedInput,
} from "@neherlab/app-contracts";
import { zCheckConfigRequest } from "@neherlab/app-contracts";
import { createApiClient } from "@neherlab/app-contracts/client";
import { Toaster } from "@neherlab/design/ui/toast";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import {
  createMemoryHistory,
  createRootRoute,
  createRoute,
  createRouter,
  RouterProvider,
} from "@tanstack/react-router";
import { configure, within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import pWaitFor from "p-wait-for";
import { act, useCallback, useEffect, useState } from "react";
import { createRoot, type Root } from "react-dom/client";
import { FormProvider, useForm, useWatch } from "react-hook-form";

import { ApiContext, type ApiContextValue } from "../../api/context";
import { useDraftStore } from "../../store/draft";
import { type FormConfig, fromFormConfig, toFormConfig } from "../formValues";
import { InputsSection } from "../InputsSection";
import { PathPicker } from "../PathPicker";
import { useDatasetLoader, useExampleIndex } from "../useExampleData";

const NATIVE_FETCH = { fetch, Request, Response, Headers, Blob, File, FormData, AbortController, AbortSignal };

const ORIGIN = "http://treetime.test";

const BASE = "/treetime/";

const WEIGHTS: string[] = [];

const MUGRATION_CONFIG = "Mugration of the zika dataset";

const ANCESTRAL_CONFIG = "Ancestral sequences of the zika dataset";

const INDEX: DataIndex = {
  files: [
    { path: "zika/20/ancestral.yaml", size: 32 },
    { path: "zika/20/aln.fasta.xz", size: 4 },
    { path: "zika/20/metadata.tsv", size: 8 },
    { path: "zika/20/mugration.yaml", size: 16 },
    { path: "zika/20/tree.nwk", size: 11 },
    { path: "zika/country_weights.csv", size: 12 },
  ],
  datasets: [
    {
      name: "zika/20",
      inputs: [
        { kind: "tree", path: "zika/20/tree.nwk" },
        { kind: "alignment", path: "zika/20/aln.fasta.xz" },
        { kind: "metadata", path: "zika/20/metadata.tsv" },
      ],
    },
  ],
  configs: [
    {
      path: "zika/20/mugration.yaml",
      command: "mugration",
      title: "Mugration of the zika dataset",
      folder: "zika/20",
      content: "tree: tree.nwk\nmetadata: metadata.tsv\nweights: ../country_weights.csv\nattribute: country\n",
      inputs: {
        metadata: ["zika/20/metadata.tsv"],
        tree: ["zika/20/tree.nwk"],
        weights: ["zika/country_weights.csv"],
      },
    },
    {
      path: "zika/20/ancestral.yaml",
      command: "ancestral",
      title: "Ancestral sequences of the zika dataset",
      folder: "zika/20",
      content: "tree: tree.nwk\nalignment: aln.fasta.xz\n",
      inputs: {
        alignment: ["zika/20/aln.fasta.xz"],
        tree: ["zika/20/tree.nwk"],
      },
    },
  ],
  pipelines: [],
};

const FILES = new Map([
  ["zika/20/aln.fasta.xz", "xz!!"],
  ["zika/20/metadata.tsv", "strain\tx\n"],
  ["zika/20/tree.nwk", "(A:1,B:1);\n"],
  ["zika/country_weights.csv", "country,w\nA,1\n"],
]);

const UPLOAD_RUN: RunRecord = {
  id: "r1",
  title: "Uploads",
  command: "timetree",
  config: {},
  status: "created",
  pinned: false,
  created_at: "2026-10-22T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [],
  changed_settings: [],
  headline: {},
  output_files: [],
  warnings: [],
};

const CODE = {
  command_line: [],
  command_line_text: "",
  yaml: [],
  yaml_text: "",
  config_file: "",
  config_command: "",
};

const mounted: Root[] = [];

/* HACK: the DOM stays registered for the rest of this file's worker (`bun test --parallel`), because React's
   scheduler runs a task queued before the last unmount that reads `window`, after an `unregister()` would remove it. */
beforeAll(() => {
  GlobalRegistrator.register();
  Object.assign(globalThis, NATIVE_FETCH);
  process.env["BASE_URL"] = BASE;
  configure({ eventWrapper: inAct });
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
  Object.assign(globalThis, { fetch: NATIVE_FETCH.fetch });
});

describe("example data of the web apps", () => {
  test("a dataset uploads its input files from the data folder under the base path", async () => {
    const server = fakeServer();
    await render(server, "timetree");
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Example datasets" }));
    await user.click(await page().findByRole("button", { name: "Load" }));
    await pWaitFor(() => Object.keys(useDraftStore.getState().draft.sources).length === 3);

    expect([
      server.static.toSorted(),
      server.uploads.toSorted(),
      settingText("Tree setting"),
      settingText("Alignment setting"),
      useDraftStore.getState().draft.sources["metadata"],
    ]).toStrictEqual([
      [
        "/treetime/data/index.json",
        "/treetime/data/zika/20/aln.fasta.xz",
        "/treetime/data/zika/20/metadata.tsv",
        "/treetime/data/zika/20/tree.nwk",
      ],
      ["aln.fasta.xz:4", "metadata.tsv:9", "tree.nwk:11"],
      JSON.stringify("/runs/r1/inputs/tree.nwk"),
      JSON.stringify(["/runs/r1/inputs/aln.fasta.xz"]),
      { label: "metadata.tsv", size: 9, origin: "upload", run_id: "r1" },
    ]);
  });

  test("a dataset loads outside the form, as the command palette loads it", async () => {
    const server = fakeServer();
    await render(server, "timetree");

    await userEvent.setup({ document }).click(page().getByRole("button", { name: "Load a dataset outside the form" }));
    await pWaitFor(() => Object.keys(useDraftStore.getState().draft.sources).length === 3);

    expect([server.uploads, settingText("Tree setting"), useDraftStore.getState().draft.sources["tree"]]).toStrictEqual(
      [
        ["tree.nwk:11", "metadata.tsv:9", "aln.fasta.xz:4"],
        JSON.stringify("/runs/r1/inputs/tree.nwk"),
        { label: "tree.nwk", size: 11, origin: "upload", run_id: "r1" },
      ],
    );
  });

  test("a dataset file that cannot be read uploads nothing", async () => {
    const server = fakeServer({ missing: "zika/20/metadata.tsv" });
    await render(server, "timetree");
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Example datasets" }));
    await user.click(await page().findByRole("button", { name: "Load" }));
    await pWaitFor(() => server.static.length === 4);

    expect([
      await page().findByRole("button", { name: "Load" }),
      server.uploads,
      settingText("Tree setting"),
    ]).toStrictEqual([page().getByRole("button", { name: "Load" }), [], "null"]);
  });

  test("a config uploads each file it reads and checks its text with the uploads as inputs", async () => {
    const server = fakeServer();
    await render(server, "timetree");

    await loadExampleConfig();

    expect([server.replaced.at(-1), useDraftStore.getState().draft.sources]).toStrictEqual([
      {
        metadata: "/runs/r1/inputs/metadata.tsv",
        tree: "/runs/r1/inputs/tree.nwk",
        weights: "/runs/r1/inputs/country_weights.csv",
      },
      {
        tree: { label: "tree.nwk", size: 11, origin: "upload", run_id: "r1" },
        metadata: { label: "metadata.tsv", size: 9, origin: "upload", run_id: "r1" },
        weights: { label: "country_weights.csv", size: 14, origin: "upload", run_id: "r1" },
      },
    ]);
  });

  test("replacing an input of a loaded example config releases its upload", async () => {
    const server = fakeServer();
    await render(server, "timetree");
    await loadExampleConfig();
    await page().findByText(JSON.stringify("/runs/r1/inputs/tree.nwk"));

    await userEvent
      .setup({ document })
      .upload(page().getByLabelText("Tree file to upload"), new File(["(C:1,D:1);\n"], "mine.nwk"));
    await pWaitFor(() => server.released.length > 0);

    expect([server.released, settingText("Tree setting")]).toStrictEqual([
      ["/runs/r1/inputs/tree.nwk"],
      JSON.stringify("/runs/r1/inputs/mine.nwk"),
    ]);
  });

  test("a file chosen while a dataset is fetched stays chosen and the dataset is dropped", async () => {
    const alignment = Promise.withResolvers<void>();
    const server = fakeServer({ held: new Map([["zika/20/aln.fasta.xz", alignment.promise]]) });
    await render(server, "timetree");
    await chooseExampleDataset();
    await pWaitFor(() => server.static.includes("/treetime/data/zika/20/aln.fasta.xz"));

    await userEvent
      .setup({ document })
      .upload(page().getByLabelText("Tree file to upload"), new File(["(C:1,D:1);\n"], "mine.nwk"));
    await page().findByText(JSON.stringify("/runs/r1/inputs/mine.nwk"));
    alignment.resolve();
    await pWaitFor(() => page().queryByRole("cell", { name: "zika/20" }) === null);

    expect([settingText("Tree setting"), settingText("Alignment setting"), server.uploads]).toStrictEqual([
      JSON.stringify("/runs/r1/inputs/mine.nwk"),
      "[]",
      ["mine.nwk:11"],
    ]);
  });

  test("an upload of an input a dataset does not set is kept when it completes after the dataset load", async () => {
    const weights = Promise.withResolvers<void>();
    const server = fakeServer({ heldUploads: new Map([["mine.csv", weights.promise]]) });
    await render(server, "mugration");

    await userEvent
      .setup({ document })
      .upload(page().getByLabelText("Weights file to upload"), new File(["country,w\nA,1\n"], "mine.csv"));
    await pWaitFor(() => server.started.includes("mine.csv"));
    await chooseExampleDataset();
    await pWaitFor(() => Object.keys(useDraftStore.getState().draft.sources).length === 2);
    weights.resolve();
    await pWaitFor(() => useDraftStore.getState().draft.sources["weights"] !== undefined);

    expect([
      settingText("Weights setting"),
      settingText("Tree setting"),
      Object.keys(useDraftStore.getState().draft.sources).toSorted(),
      server.released,
    ]).toStrictEqual([
      JSON.stringify("/runs/r1/inputs/mine.csv"),
      JSON.stringify("/runs/r1/inputs/tree.nwk"),
      ["metadata", "tree", "weights"],
      [],
    ]);
  });

  test("a dataset chosen while an example config is fetched replaces the config and releases its uploads", async () => {
    const weights = Promise.withResolvers<void>();
    const server = fakeServer({ held: new Map([["zika/country_weights.csv", weights.promise]]) });
    await render(server, "timetree");
    await chooseExampleConfig(MUGRATION_CONFIG);
    await pWaitFor(() => server.static.includes("/treetime/data/zika/country_weights.csv"));

    await chooseExampleDataset();
    await pWaitFor(() => Object.keys(useDraftStore.getState().draft.sources).length === 3);
    weights.resolve();
    await pWaitFor(() => server.released.length === 2);

    expect([
      useDraftStore.getState().draft.command,
      server.released.toSorted(),
      settingText("Tree setting"),
    ]).toStrictEqual([
      "timetree",
      ["/runs/r1/inputs/metadata.tsv", "/runs/r1/inputs/tree.nwk"],
      JSON.stringify("/runs/r1/inputs/again/tree.nwk"),
    ]);
  });

  test("a config of another command releases every upload of the earlier draft that it does not keep", async () => {
    const server = fakeServer();
    await render(server, "timetree");
    await loadExampleConfig();

    await loadExampleConfig(ANCESTRAL_CONFIG, "ancestral");
    await pWaitFor(() => server.released.length === 3);

    expect(server.released.toSorted()).toStrictEqual([
      "/runs/r1/inputs/country_weights.csv",
      "/runs/r1/inputs/metadata.tsv",
      "/runs/r1/inputs/tree.nwk",
    ]);
  });

  test("a config uploads its inputs as replacements of the uploads of the draft it replaces", async () => {
    const server = fakeServer();
    await render(server, "timetree");
    await loadExampleConfig();

    await loadExampleConfig(ANCESTRAL_CONFIG, "ancestral");

    expect(server.replaces.slice(3).map((paths) => paths.toSorted())).toStrictEqual(
      Array.from({ length: 2 }, () => [
        "/runs/r1/inputs/country_weights.csv",
        "/runs/r1/inputs/metadata.tsv",
        "/runs/r1/inputs/tree.nwk",
      ]),
    );
  });

  test("a config the check rejects releases its uploads and keeps the form", async () => {
    const server = fakeServer({ rejected: true });
    await render(server, "timetree");
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Example configs" }));
    await user.click(await loadButton(MUGRATION_CONFIG));
    await pWaitFor(() => server.released.length === 3);

    expect([server.released.toSorted(), useDraftStore.getState().draft]).toStrictEqual([
      ["/runs/r1/inputs/country_weights.csv", "/runs/r1/inputs/metadata.tsv", "/runs/r1/inputs/tree.nwk"],
      expect.objectContaining({ command: "timetree", sources: {} }),
    ]);
  });
});

async function loadExampleConfig(title = MUGRATION_CONFIG, command: AppCommand = "mugration") {
  await chooseExampleConfig(title);
  await pWaitFor(() => useDraftStore.getState().draft.command === command);
}

async function chooseExampleConfig(title: string) {
  const user = userEvent.setup({ document });

  if (page().queryByRole("button", { name: "Example configs" }) !== null) {
    await user.click(page().getByRole("button", { name: "Example configs" }));
  }

  await user.click(await loadButton(title));
}

async function chooseExampleDataset() {
  const user = userEvent.setup({ document });
  await user.click(page().getByRole("button", { name: "Example datasets" }));
  await user.click(await page().findByRole("button", { name: "Load" }));
}

async function loadButton(title: string): Promise<HTMLElement> {
  return within(await page().findByRole("row", { name: new RegExp(title, "u") })).getByRole("button", { name: "Load" });
}

function fakeServer({
  missing,
  rejected = false,
  held = new Map(),
  heldUploads = new Map(),
}: {
  missing?: string;
  rejected?: boolean;
  held?: ReadonlyMap<string, Promise<void>>;
  heldUploads?: ReadonlyMap<string, Promise<void>>;
} = {}) {
  const requested: string[] = [];
  const uploads: string[] = [];
  const released: string[] = [];
  const replaced: SparseConfig[] = [];
  const replaces: string[][] = [];
  const stored = new Set<string>();
  const started: string[] = [];

  const fetch = async (input: Request | URL | string, init?: RequestInit): Promise<Response> => {
    const request = new Request(input instanceof Request ? input : new URL(input, ORIGIN), init);
    const url = new URL(request.url);
    const { pathname } = url;
    const file = pathname.startsWith(`${BASE}data/`) ? decodeURIComponent(pathname.slice(BASE.length + 5)) : undefined;

    if (file !== undefined) {
      requested.push(pathname);
      const text = FILES.get(file);

      if (file === "index.json") {
        return Response.json(INDEX);
      }

      await held.get(file);

      return text === undefined || file === missing ? new Response("not found", { status: 404 }) : new Response(text);
    }

    if (request.method === "GET" && pathname === "/api/runs") {
      return Response.json({ runs: [], active_runs: 0 });
    }

    if (request.method === "POST" && pathname === "/api/runs") {
      return Response.json(UPLOAD_RUN);
    }

    if (request.method === "GET" && pathname === "/api/runs/r1") {
      return Response.json(UPLOAD_RUN);
    }

    if (request.method === "DELETE" && pathname === "/api/runs/r1/inputs") {
      released.push(url.searchParams.get("path") ?? "");

      return new Response(null, { status: 200 });
    }

    const name = /^\/api\/runs\/r1\/inputs\/(?<name>[^/]+)$/u.exec(pathname)?.groups?.["name"];

    if (request.method === "PUT" && name !== undefined) {
      started.push(name);
      await heldUploads.get(name);
      const size = (await request.arrayBuffer()).byteLength;
      uploads.push(`${name}:${size}`);
      replaces.push(url.searchParams.getAll("replaces"));
      const path = stored.has(name) ? `/runs/r1/inputs/again/${name}` : `/runs/r1/inputs/${name}`;
      stored.add(name);
      const uploaded: UploadedInput = { name, path, size, sha256: "0" };

      return Response.json(uploaded);
    }

    if (request.method === "POST" && pathname === "/api/check-config") {
      const body = zCheckConfigRequest.parse(await request.json());
      replaced.push(body.replace_inputs ?? {});

      return Response.json(checkAnswer(body.command, body.replace_inputs ?? {}, rejected));
    }

    return Response.json({ message: `unexpected ${request.method} ${pathname}` }, { status: 500 });
  };

  return { fetch, static: requested, uploads, released, replaced, replaces, started };
}

function checkAnswer(command: AppCommand, replaceInputs: SparseConfig, rejected: boolean): CheckConfigResponse {
  if (rejected) {
    return {
      status: "invalid",
      command,
      message: "the weights do not parse",
      causes: [],
      problems: [],
      messages: ["the weights do not parse"],
      checks: [],
      choices: [],
    };
  }

  return {
    status: "valid",
    command,
    config: command === "mugration" ? { ...replaceInputs, attribute: "country" } : replaceInputs,
    code: CODE,
    checks: [],
    choices: [],
  };
}

async function render(server: ReturnType<typeof fakeServer>, command: AppCommand) {
  Object.assign(globalThis, { fetch: server.fetch });
  useDraftStore.getState().reset(command);
  const api = apiContext(server.fetch);
  const queryClient = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  const root = createRootRoute({ component: DraftPage });

  const router = createRouter({
    routeTree: root.addChildren([
      createRoute({ getParentRoute: () => root, path: "/" }),
      createRoute({ getParentRoute: () => root, path: "/new" }),
    ]),
    history: createMemoryHistory({ initialEntries: ["/"] }),
  });

  const container = document.createElement("div");
  document.body.append(container);
  const reactRoot = createRoot(container);
  mounted.push(reactRoot);
  inAct(() =>
    reactRoot.render(
      <ApiContext.Provider value={api}>
        <QueryClientProvider client={queryClient}>
          <Toaster>
            <RouterProvider router={router} />
          </Toaster>
        </QueryClientProvider>
      </ApiContext.Provider>,
    ),
  );
  await page().findByRole("button", { name: "Example datasets" });
}

function apiContext(fetch: typeof globalThis.fetch): ApiContextValue {
  return { client: createApiClient({ baseUrl: ORIGIN, fetch }) };
}

function DraftPage() {
  const epoch = useDraftStore((state) => state.epoch);
  const command = useDraftStore((state) => state.draft.command);

  return (
    <>
      <DatasetOutsideForm command={command} />
      <DraftForm key={`${command}:${epoch}`} command={command} />
    </>
  );
}

function DatasetOutsideForm({ command }: { command: AppCommand }) {
  const { data: examples } = useExampleIndex();
  const loadDataset = useDatasetLoader();
  const dataset = examples?.index?.datasets[0];

  const onClick = useCallback(() => {
    if (examples !== undefined && dataset !== undefined) {
      void loadDataset(examples, dataset, command);
    }
  }, [command, dataset, examples, loadDataset]);

  return (
    <button type="button" onClick={onClick}>
      Load a dataset outside the form
    </button>
  );
}

function DraftForm({ command }: { command: AppCommand }) {
  const [initial] = useState(() => toFormConfig(useDraftStore.getState().draft.config));
  const form = useForm<FormConfig>({ defaultValues: initial });

  useEffect(
    () =>
      form.subscribe({
        formState: { values: true },
        callback: ({ values }) => useDraftStore.getState().setConfig(fromFormConfig(values)),
      }),
    [form],
  );

  return (
    <FormProvider {...form}>
      <InputsSection command={command} facts={undefined} />
      <PathPicker
        command={command}
        settingKey="weights"
        title="Weights"
        extensions={WEIGHTS}
        list={false}
        filled={false}
      />
      <SettingValue name="tree" label="Tree setting" />
      <SettingValue name="alignment" label="Alignment setting" />
      <SettingValue name="weights" label="Weights setting" />
    </FormProvider>
  );
}

function SettingValue({ name, label }: { name: string; label: string }) {
  const value = useWatch<FormConfig>({ name });

  return <output aria-label={label}>{JSON.stringify(value ?? null)}</output>;
}

function inAct(callback: () => void) {
  Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

  try {
    act(callback);
  } finally {
    Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: false });
  }
}

function settingText(label: string): string | null {
  return page().getByLabelText(label).textContent;
}

function page() {
  return within(document.body);
}
