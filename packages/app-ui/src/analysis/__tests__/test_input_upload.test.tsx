import { GlobalRegistrator } from "@happy-dom/global-registrator";
import type { RunRecord, UploadedInput } from "@neherlab/app-contracts";
import { createApiClient } from "@neherlab/app-contracts/client";
import { configure, fireEvent, within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import pWaitFor from "p-wait-for";
import { act } from "react";
import { createRoot, type Root } from "react-dom/client";
import { FormProvider, useForm, useWatch } from "react-hook-form";

import { ApiContext, type ApiContextValue } from "../../api/context";
import { useDraftStore } from "../../store/draft";
import { Toaster } from "../../ui/toast";
import type { FormConfig } from "../formValues";
import { InputsSection } from "../InputsSection";

const NATIVE_FETCH = { fetch, Request, Response, Headers, Blob, File, FormData, AbortController, AbortSignal };

const UPLOAD_RUN: RunRecord = {
  id: "r1",
  title: "Uploads",
  command: "ancestral",
  config: {},
  status: "created",
  pinned: false,
  created_at: "2026-10-18T08:00:00Z",
  treetime_version: "1.0.0",
  inputs: [],
  changed_settings: [],
  headline: {},
  output_files: [],
  warnings: [],
};

const mounted: Root[] = [];

/* HACK: the DOM stays registered for the rest of this file's worker (`bun test --parallel`), because React's
   scheduler runs a task queued before the last unmount that reads `window`, after an `unregister()` would remove it. */
beforeAll(() => {
  GlobalRegistrator.register();
  Object.assign(globalThis, NATIVE_FETCH);
  configure({ eventWrapper: inAct });
});

beforeEach(() => {
  useDraftStore.getState().reset("ancestral");
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("alignment inputs of the web app", () => {
  test("uploading two files with the same name stores both in selection order", async () => {
    const server = fakeServer();
    await render(server.fetch);

    await userEvent
      .setup({ document })
      .upload(page().getByLabelText("Alignment files to upload"), [
        fasta("aligned.fasta", ">A\nACGT\n"),
        fasta("aligned.fasta", ">A\nAC\n"),
      ]);

    expect(await settingShowing(["/runs/r1/inputs/1/aligned.fasta", "/runs/r1/inputs/2/aligned.fasta"])).toBe(
      page().getByLabelText("Alignment setting"),
    );
    expect(server.uploads).toStrictEqual(["aligned.fasta:8", "aligned.fasta:6"]);
  });

  test("a new selection replaces the list", async () => {
    const server = fakeServer();
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");

    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n"), fasta("na.fasta", ">A\nAC\n")]);
    expect(await settingShowing(["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"])).toBe(
      page().getByLabelText("Alignment setting"),
    );
    await user.upload(input, [fasta("mp.fasta", ">A\nA\n")]);

    expect(await settingShowing(["/runs/r1/inputs/3/mp.fasta"])).toBe(page().getByLabelText("Alignment setting"));
  });

  test("replacing two uploaded files releases both earlier uploads on the server", async () => {
    const server = fakeServer();
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n"), fasta("na.fasta", ">A\nAC\n")]);
    await settingShowing(["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"]);

    await user.upload(input, [fasta("mp.fasta", ">A\nA\n")]);
    await pWaitFor(() => server.released.length === 2);

    expect(server.released.toSorted()).toStrictEqual(["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"]);
  });

  test("an upload that fails part-way releases the files it already stored", async () => {
    const server = fakeServer({ failing: "na.fasta" });
    await render(server.fetch);

    await userEvent
      .setup({ document })
      .upload(page().getByLabelText("Alignment files to upload"), [
        fasta("ha.fasta", ">A\nACGT\n"),
        fasta("na.fasta", ">A\nAC\n"),
      ]);
    await pWaitFor(() => server.released.length === 1);

    expect([server.released, page().getByLabelText("Alignment setting").textContent]).toStrictEqual([
      ["/runs/r1/inputs/1/ha.fasta"],
      "null",
    ]);
  });

  test("replacing files of a started run sends no release and shows no error", async () => {
    const server = fakeServer();
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n")]);
    await settingShowing(["/runs/r1/inputs/1/ha.fasta"]);
    inAct(() => useDraftStore.getState().clearRuns());

    await user.upload(input, [fasta("mp.fasta", ">A\nA\n")]);
    await settingShowing(["/runs/r2/inputs/2/mp.fasta"]);

    expect([server.released, page().queryByText("An earlier upload cannot be removed")]).toStrictEqual([[], null]);
  });

  test("removing the files during an upload keeps the setting empty", async () => {
    const server = fakeServer({ held: "na.fasta" });
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n")]);
    await settingShowing(["/runs/r1/inputs/1/ha.fasta"]);
    await user.upload(input, [fasta("na.fasta", ">A\nAC\n")]);
    await pWaitFor(() => server.holding());

    await user.click(page().getByRole("button", { name: "Remove" }));
    server.resume();
    await pWaitFor(() => server.released.includes("/runs/r1/inputs/2/na.fasta") || settingText() !== "[]");

    expect([settingText(), server.released.toSorted()]).toStrictEqual([
      "[]",
      ["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"],
    ]);
  });

  test("of two overlapping uploads the newer selection wins", async () => {
    const server = fakeServer({ held: "ha.fasta" });
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n")]);
    await pWaitFor(() => server.holding());

    await user.upload(input, [fasta("na.fasta", ">A\nAC\n")]);
    await settingShowing(["/runs/r1/inputs/1/na.fasta"]);
    server.resume();
    await pWaitFor(
      () => server.released.length > 0 || settingText() !== JSON.stringify(["/runs/r1/inputs/1/na.fasta"]),
    );

    expect([settingText(), server.released]).toStrictEqual([
      JSON.stringify(["/runs/r1/inputs/1/na.fasta"]),
      ["/runs/r1/inputs/2/ha.fasta"],
    ]);
  });

  test("an upload answer that arrives after a start leaves the setting unchanged", async () => {
    const server = fakeServer({ held: "na.fasta" });
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n")]);
    await settingShowing(["/runs/r1/inputs/1/ha.fasta"]);
    await user.upload(input, [fasta("na.fasta", ">A\nAC\n")]);
    await pWaitFor(() => server.holding());

    inAct(() => useDraftStore.getState().clearRuns());
    server.resume();
    await page().findByRole("button", { name: "Replace" });

    expect([settingText(), server.released, page().queryByText("An earlier upload cannot be removed")]).toStrictEqual([
      JSON.stringify(["/runs/r1/inputs/1/ha.fasta"]),
      [],
      null,
    ]);
  });

  test.each([
    ["a reset", () => useDraftStore.getState().reset("ancestral")],
    ["a start", () => useDraftStore.getState().clearRuns()],
  ])("an upload run created after %s receives no files and assigns nothing", async (_, discard) => {
    const server = fakeServer({ holdCreates: true });
    await render(server.fetch);
    const user = userEvent.setup({ document });
    await user.upload(page().getByLabelText("Alignment files to upload"), [fasta("ha.fasta", ">A\nACGT\n")]);
    await pWaitFor(() => server.creates() > 0);

    inAct(discard);
    server.resumeCreatesInReverse();
    await page().findByRole("button", { name: "Choose files" });

    expect([
      useDraftStore.getState().draft.upload_run_id,
      settingText(),
      server.uploads,
      server.released,
    ]).toStrictEqual([undefined, "null", [], []]);
  });

  test("two first uploads at the same time share one upload run and both end assigned", async () => {
    const server = fakeServer({ holdCreates: true });
    await render(server.fetch);
    const user = userEvent.setup({ document });

    await user.upload(page().getByLabelText("Tree file to upload"), [fasta("tree.nwk", "(A,B);\n")]);
    await user.upload(page().getByLabelText("Alignment files to upload"), [fasta("ha.fasta", ">A\nACGT\n")]);
    await pWaitFor(() => server.creates() > 0);
    server.resumeCreatesInReverse();

    expect(await page().findByText(JSON.stringify(["/runs/r1/inputs/2/ha.fasta"]))).toBe(
      page().getByLabelText("Alignment setting"),
    );
    expect([page().getByLabelText("Tree setting").textContent, server.creates()]).toStrictEqual([
      JSON.stringify("/runs/r1/inputs/1/tree.nwk"),
      1,
    ]);
  });

  test("an upload that replaces uploaded files names them", async () => {
    const server = fakeServer();
    await render(server.fetch);
    const user = userEvent.setup({ document });
    const input = page().getByLabelText("Alignment files to upload");
    await user.upload(input, [fasta("ha.fasta", ">A\nACGT\n"), fasta("na.fasta", ">A\nAC\n")]);
    await settingShowing(["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"]);

    await user.upload(input, [fasta("mp.fasta", ">A\nA\n")]);
    await settingShowing(["/runs/r1/inputs/3/mp.fasta"]);

    expect(server.replaces).toStrictEqual([[], [], ["/runs/r1/inputs/1/ha.fasta", "/runs/r1/inputs/2/na.fasta"]]);
  });

  test("dropping two files on the alignment row stores both in drop order", async () => {
    const server = fakeServer();
    await render(server.fetch);
    const files = [fasta("na.fasta", ">A\nAC\n"), fasta("ha.fasta", ">A\nACGT\n")];

    fireEvent.drop(page().getByRole("group", { name: "Alignment" }), {
      dataTransfer: {
        files,
        items: files.map((file) => ({ kind: "file", type: file.type, getAsFile: () => file })),
        types: ["Files"],
      },
    });

    expect(await settingShowing(["/runs/r1/inputs/1/na.fasta", "/runs/r1/inputs/2/ha.fasta"])).toBe(
      page().getByLabelText("Alignment setting"),
    );
  });
});

function fakeServer({ failing, held, holdCreates }: { failing?: string; held?: string; holdCreates?: boolean } = {}) {
  const uploads: string[] = [];
  const released: string[] = [];
  const replaces: string[][] = [];
  const gate = Promise.withResolvers<undefined>();
  let holding = false;

  let runs = 0;
  const heldCreates: (() => void)[] = [];
  let createsOpen = holdCreates !== true;

  const fetch = async (input: Request | URL | string, init?: RequestInit): Promise<Response> => {
    const request = new Request(input, init);
    const { pathname } = new URL(request.url);

    if (request.method === "POST" && pathname === "/api/runs") {
      runs += 1;
      const record = { ...UPLOAD_RUN, id: `r${runs}` };

      if (!createsOpen) {
        const answer = Promise.withResolvers<undefined>();
        heldCreates.push(() => answer.resolve(undefined));
        await answer.promise;
      }

      return Response.json(record);
    }

    const run = /^\/api\/runs\/(?<id>r\d+)/u.exec(pathname)?.groups?.["id"] ?? "";

    if (request.method === "GET" && pathname === `/api/runs/${run}`) {
      return Response.json({ ...UPLOAD_RUN, id: run });
    }

    const name = /^\/api\/runs\/r\d+\/inputs\/(?<name>[^/]+)$/u.exec(pathname)?.groups?.["name"];

    if (request.method === "DELETE" && pathname === `/api/runs/${run}/inputs`) {
      released.push(new URL(request.url).searchParams.get("path") ?? "");

      return run === `r${runs}`
        ? new Response(null, { status: 200 })
        : Response.json({ code: "conflict", message: "the run has started", causes: [] }, { status: 409 });
    }

    if (request.method === "PUT" && name === failing) {
      return Response.json({ code: "too-large", message: "the inputs do not fit", causes: [] }, { status: 413 });
    }

    if (request.method === "PUT" && name !== undefined) {
      replaces.push(new URL(request.url).searchParams.getAll("replaces"));

      if (name === held) {
        holding = true;
        await gate.promise;
      }

      const size = (await request.arrayBuffer()).byteLength;

      uploads.push(`${name}:${size}`);

      const uploaded: UploadedInput = {
        name,
        path: `/runs/${run}/inputs/${uploads.length}/${name}`,
        size,
        sha256: "0",
      };

      return Response.json(uploaded);
    }

    return Response.json({ message: `unexpected ${request.method} ${pathname}` }, { status: 500 });
  };

  const resumeCreatesInReverse = () => {
    createsOpen = true;

    for (const answer of heldCreates.toReversed()) {
      answer();
    }
  };

  return {
    uploads,
    released,
    replaces,
    fetch,
    holding: () => holding,
    resume: () => gate.resolve(undefined),
    creates: () => runs,
    resumeCreatesInReverse,
  };
}

async function render(fetch: typeof globalThis.fetch) {
  const api = apiContext(fetch);
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() =>
    root.render(
      <ApiContext.Provider value={api}>
        <Toaster>
          <Form />
        </Toaster>
      </ApiContext.Provider>,
    ),
  );
  await page().findByRole("group", { name: "Alignment" });
}

function apiContext(fetch: typeof globalThis.fetch): ApiContextValue {
  return { client: createApiClient({ baseUrl: "http://localhost", fetch }) };
}

function inAct(callback: () => void) {
  Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

  try {
    act(callback);
  } finally {
    Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: false });
  }
}

async function settingShowing(paths: string[]): Promise<HTMLElement> {
  return page().findByText(JSON.stringify(paths));
}

function settingText(): string | null {
  return page().getByLabelText("Alignment setting").textContent;
}

function page() {
  return within(document.body);
}

function fasta(name: string, text: string): File {
  return new NATIVE_FETCH.File([text], name, { type: "text/plain" });
}

function Form() {
  const form = useForm<FormConfig>();

  return (
    <FormProvider {...form}>
      <InputsSection command="ancestral" facts={undefined} />
      <SettingValue name="tree" label="Tree setting" />
      <SettingValue name="alignment" label="Alignment setting" />
    </FormProvider>
  );
}

function SettingValue({ name, label }: { name: string; label: string }) {
  const value = useWatch<FormConfig>({ name });

  return <output aria-label={label}>{JSON.stringify(value ?? null)}</output>;
}
