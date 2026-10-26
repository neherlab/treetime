import "../../__tests__/register-dom";
import type { CommandSettings, SettingSpec } from "@neherlab/app-contracts";
import { LocaleProvider } from "@neherlab/app-ui/i18n";
import { userEvent } from "@testing-library/user-event";
import pWaitFor from "p-wait-for";
import { act } from "react";

import { page, render, setupDom } from "../../__tests__/render";
import type { SearchPages, SearchResult } from "../../lib/pagefind";
import { manifestOf, siteVersions } from "../../lib/version-state";
import { DocsSearch } from "../DocsSearch";
import { SettingsReference } from "../SettingsReference";
import { VersionBanner } from "../VersionBanner";

setupDom();

describe("docs search", () => {
  test("opens with Mod+K and lists titles and highlighted excerpts of the results", async () => {
    render(
      <LocaleProvider locale="en">
        <DocsSearch loadSearch={loadFakeSearch} />
      </LocaleProvider>,
    );
    const user = userEvent.setup({ document });

    await user.keyboard("{Control>}k{/Control}");
    await user.type(page().getByRole("combobox", { name: "Search the docs" }), "clock");
    await pWaitFor(() => page().queryByText("Math") !== null);
    const marked = [...document.querySelectorAll("mark")].map((mark) => mark.textContent);

    expect([page().getByText("Math").textContent, marked]).toStrictEqual(["Math", ["clock"]]);
  });

  test("leaves Shift+Home and Shift+End in the search box to the browser, which extends the text selection", async () => {
    render(
      <LocaleProvider locale="en">
        <DocsSearch loadSearch={loadFakeSearch} />
      </LocaleProvider>,
    );
    const user = userEvent.setup({ document });
    const keys: KeyboardEvent[] = [];
    const recording = new AbortController();

    await user.keyboard("{Control>}k{/Control}");
    await user.type(page().getByRole("combobox", { name: "Search the docs" }), "clock");
    await pWaitFor(() => page().queryByText("Math") !== null);
    window.addEventListener(
      "keydown",
      (event) => {
        if (event.key !== "Shift") {
          keys.push(event);
        }
      },
      { capture: true, signal: recording.signal },
    );
    await user.keyboard("{Shift>}{Home}{/Shift}{Shift>}{End}{/Shift}");
    recording.abort();

    expect(keys.map((event) => [event.key, event.defaultPrevented])).toStrictEqual([
      ["Home", false],
      ["End", false],
    ]);
  });

  test("says that search works on the built site when the index does not load", async () => {
    render(
      <LocaleProvider locale="en">
        <DocsSearch loadSearch={loadMissingSearch} />
      </LocaleProvider>,
    );
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Search the docs" }));
    await pWaitFor(() => page().queryByText(/Search works on the built site/u) !== null);

    expect(page().getByText(/Search works on the built site/u).tagName).toBe("P");
  });
});

describe("settings reference", () => {
  test("lists every setting under its group heading, in the order of the help", () => {
    render(<SettingsReference settings={SETTINGS} locale="en" />);

    expect(
      page()
        .getAllByRole("heading", { level: 2 })
        .map((heading) => heading.textContent),
    ).toStrictEqual(["Inputs", "Input data", "Molecular clock"]);
    expect(
      page()
        .getAllByRole("heading", { level: 3 })
        .map((heading) => heading.textContent),
    ).toStrictEqual(["Tree", "Clock rate", "Clock filter"]);
  });

  test("filters the settings by flag and announces the count, and clearing the filter restores them", async () => {
    render(<SettingsReference settings={SETTINGS} locale="en" />);
    const user = userEvent.setup({ document });
    const filter = page().getByRole("searchbox", { name: "Filter settings" });

    await user.type(filter, "--clock-rate");
    const filtered = [entryKeys(), page().getByRole("status").textContent];
    await user.type(filter, "{Control>}a{/Control}{Backspace}");

    expect([filtered, [entryKeys(), page().getByRole("status").textContent]]).toStrictEqual([
      [["clock_rate"], "1 setting matches"],
      [["tree", "clock_rate", "clock_filter"], "3 settings match"],
    ]);
  });
});

describe("version banner", () => {
  test.each([
    ["shows nothing on the preferred version", "dev", null],
    [
      "names the preferred version on an older one",
      "1.0",
      "These are the docs of TreeTime 1.0. The recommended version is dev.",
    ],
  ])("%s", async (_, version, expected) => {
    installFetch({
      "/versions.json": [
        { version: "dev", label: "dev", preferred: true },
        { version: "1.0", label: "1.0" },
      ],
      "/dev/pages.json": { version: "dev", pages: { en: ["", "showcase/math"] } },
    });
    render(<VersionBanner version={version} page="showcase/math" locale="en" />);
    await act(async () => {
      await siteVersions();
      await manifestOf("dev");
    });
    const banner = page().queryByRole("note");

    expect([
      banner?.textContent?.split(" Open")[0] ?? null,
      banner?.querySelector("a")?.getAttribute("href") ?? null,
    ]).toStrictEqual(expected === null ? [null, null] : [expected, "/dev/en/showcase/math"]);
  });
});

const RESULTS: SearchResult[] = [
  {
    url: "/dev/en/showcase/math",
    excerpt: "The <mark>clock</mark> rate &amp; dates",
    meta: { title: "Math" },
    sub_results: [],
  },
];

function loadFakeSearch(): Promise<SearchPages> {
  return Promise.resolve(searchFake);
}

function searchFake(): Promise<SearchResult[]> {
  return Promise.resolve(RESULTS);
}

function loadMissingSearch(): Promise<SearchPages> {
  return Promise.reject(new Error("pagefind.js is missing"));
}

function installFetch(responses: Record<string, unknown>): void {
  Object.assign(globalThis, {
    fetch: (url: string) => {
      const body = responses[url];

      return Promise.resolve(body === undefined ? new Response(null, { status: 404 }) : Response.json(body));
    },
  });
}

function entryKeys(): string[] {
  return [...document.querySelectorAll("article")].map((entry) => entry.id);
}

const TREE: SettingSpec = spec({ key: "tree", label: "Tree", flag: "--tree", group: "Input data", role: "input" });

const CLOCK_RATE: SettingSpec = spec({
  key: "clock_rate",
  label: "Clock rate",
  flag: "--clock-rate",
  group: "Molecular clock",
});

const CLOCK_FILTER: SettingSpec = spec({
  key: "clock_filter",
  label: "Clock filter",
  flag: "--clock-filter",
  group: "Molecular clock",
});

const SETTINGS: CommandSettings = {
  command: "clock",
  title: "Clock signal",
  description: "Estimate the clock rate.",
  inputs: [
    { kind: "tree", need: "required", label: "Tree", formats: "Newick or Nexus", extensions: ["nwk"], list: false },
  ],
  uses_dates: true,
  groups: ["Input data", "Molecular clock"],
  main_settings: [],
  choices: [],
  settings: [TREE, CLOCK_RATE, CLOCK_FILTER],
};

function spec(fields: Pick<SettingSpec, "key" | "label" | "flag" | "group"> & Partial<SettingSpec>): SettingSpec {
  return {
    role: "setting",
    kind: "number",
    nullable: true,
    options: [],
    item_kind: "string",
    fields: [],
    examples: [],
    value_names: ["VALUE"],
    conflicts: [],
    help: `${fields.label} of the run`,
    more: "",
    ...fields,
  };
}
