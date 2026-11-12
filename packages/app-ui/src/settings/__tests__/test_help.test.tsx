import type { JsonValue, SettingSpec } from "@neherlab/app-contracts";
import { within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import type { ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import { LocaleProvider, type Locale } from "../../i18n/locale";
import { commandLineText, defaultText, SettingHelp } from "../help";

const TREE_DIALECT: SettingSpec = {
  key: "newick_dialect",
  label: "Newick dialect",
  flag: "--newick-dialect",
  group: "Input data",
  purpose: "data",
  set_by_app: false,
  role: "setting",
  kind: "text",
  nullable: false,
  options: [],
  item_kind: "string",
  fields: [],
  default_value: "classic,beast",
  examples: [],
  value_names: ["DIALECT"],
  conflicts: [],
  help: "Dialect of the tree as `<structure>,<annotations>`, for example <structure>",
  more: "Annotations such as `[&&NHX:key=value]` are read.\n\nA second paragraph.",
};

const mounted: Root[] = [];

beforeAll(registerDom);

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("setting help", () => {
  test("renders code spans as code and keeps angle-bracket words as text", () => {
    render(<SettingHelp spec={TREE_DIALECT} />, "en");

    expect([
      page().getByText("<structure>,<annotations>").tagName,
      page().getByText(/for example <structure>$/u).tagName,
    ]).toStrictEqual(["CODE", "SPAN"]);
  });

  test("holds the extra text in the page before More opens it", async () => {
    render(<SettingHelp spec={TREE_DIALECT} />, "en");
    const user = userEvent.setup({ document });
    const extra = page().getByText("A second paragraph.");
    const before = [document.body.contains(extra), extra.closest("[hidden]") !== null];

    await user.click(page().getByRole("button", { name: "More" }));

    expect([before, extra.closest("[hidden]") !== null]).toStrictEqual([[true, true], false]);
  });

  test.each<[Locale, string]>([
    ["en", "More"],
    ["de", "Mehr"],
  ])("names the extra text toggle in %s", (locale, name) => {
    render(<SettingHelp spec={TREE_DIALECT} />, locale);

    expect(page().getByRole("button").textContent).toBe(name);
  });
});

describe("command line of a setting", () => {
  test.each<[string, SettingSpec, JsonValue | undefined, string]>([
    ["a value", { ...TREE_DIALECT, kind: "enum" }, "rich,nhx", "--newick-dialect rich,nhx"],
    ["a number", { ...TREE_DIALECT, flag: "--clock-filter", kind: "number" }, 4, "--clock-filter 4"],
    ["a list", { ...TREE_DIALECT, flag: "--relax", kind: "list" }, [1, 0.5], "--relax 1 0.5"],
    ["a switch", { ...TREE_DIALECT, flag: "--keep-root", kind: "switch" }, true, "--keep-root"],
    ["no value", TREE_DIALECT, undefined, "--newick-dialect DIALECT"],
  ])("shows %s", (_name, spec, value, expected) => {
    expect(commandLineText(spec, value)).toBe(expected);
  });
});

describe("default text", () => {
  test.each<[string, JsonValue | undefined, Locale, string]>([
    ["an unset value in en", undefined, "en", "not set"],
    ["an unset value in de", undefined, "de", "nicht gesetzt"],
    ["an empty list in en", [], "en", "none"],
    ["an empty list in de", [], "de", "keine"],
    ["a list in en", ["a", "b"], "en", "a b"],
    ["true in en", true, "en", "on"],
    ["true in de", true, "de", "an"],
    ["false in en", false, "en", "off"],
    ["false in de", false, "de", "aus"],
    ["a number in de", 0.001, "de", "0.001"],
    ["an object in en", { name: "cp3" }, "en", '{"name":"cp3"}'],
  ])("writes %s", (_, value, locale, expected) => {
    expect(defaultText(value, locale)).toBe(expected);
  });
});

function render(node: ReactNode, locale: Locale) {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(<LocaleProvider locale={locale}>{node}</LocaleProvider>));
}

function page() {
  return within(document.body);
}
