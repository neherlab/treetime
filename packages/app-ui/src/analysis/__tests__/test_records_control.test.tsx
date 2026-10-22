import type { SettingSpec } from "@neherlab/app-contracts";
import { within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import { createRoot, type Root } from "react-dom/client";
import { FormProvider, useForm, useWatch } from "react-hook-form";

import { inAct, registerDom } from "../../__tests__/dom";
import type { FormConfig, FormValue } from "../formValues";
import { SettingControl } from "../SettingControl";

const ENTRIES: SettingSpec = {
  key: "entries",
  label: "Entries",
  flag: "--entries",
  group: "Input data",
  role: "setting",
  kind: "records",
  nullable: false,
  options: [],
  item_kind: "string",
  fields: [
    {
      key: "name",
      label: "Name",
      kind: "text",
      nullable: false,
      options: [],
      item_kind: "string",
      examples: [],
      help: "Name of the entry",
      more: "",
    },
    {
      key: "sites",
      label: "Sites",
      kind: "list",
      nullable: true,
      options: [],
      item_kind: "string",
      examples: [],
      help: "Sites of the entry",
      more: "",
    },
    {
      key: "model",
      label: "Model",
      kind: "enum",
      nullable: true,
      options: [
        { value: "infer", help: "" },
        { value: "jc69", help: "" },
      ],
      item_kind: "string",
      documented_default: "infer",
      examples: [],
      help: "Model of the entry",
      more: "",
    },
  ],
  examples: [],
  value_names: ["ENTRY"],
  conflicts: [],
  help: "Entries",
  more: "",
};

const mounted: Root[] = [];

beforeAll(registerDom);

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("records control", () => {
  test("shows one row per entry with the value of each field", () => {
    render([
      { name: "cp1", sites: ["1-9\\3", "20"] },
      { name: "rest", model: "jc69" },
    ]);

    expect(
      page()
        .getAllByRole("row")
        .map((row) => [
          within(row).getByRole<HTMLInputElement>("textbox", { name: /^Name/u }).value,
          within(row).getByRole<HTMLInputElement>("textbox", { name: /^Sites/u }).value,
          within(row).getByRole<HTMLSelectElement>("combobox", { name: /^Model/u }).value,
        ]),
    ).toStrictEqual([
      ["cp1", "1-9\\3 20", ""],
      ["rest", "", "jc69"],
    ]);
  });

  test("an unset field shows its documented default", () => {
    render([{ name: "cp1" }]);
    const model = page().getByRole<HTMLSelectElement>("combobox", { name: "Model of Entries entry 1" });

    expect(model.selectedOptions[0]?.textContent).toBe("Default: infer");
  });

  test("shows that the list has no entries", () => {
    render([]);

    expect(page().getByRole("grid", { name: "Entries" }).textContent).toBe("No entries");
  });

  test("editing the fields of an entry writes them to that entry", async () => {
    render([{ name: "cp1" }, { name: "rest" }]);
    const user = userEvent.setup({ document });

    await user.tripleClick(page().getByRole("textbox", { name: "Name of Entries entry 2" }));
    await user.keyboard("noncoding");
    await user.type(page().getByRole("textbox", { name: "Sites of Entries entry 2" }), "101-1600\\3 7");
    await user.selectOptions(page().getByRole("combobox", { name: "Model of Entries entry 2" }), "jc69");

    expect(settingValue()).toStrictEqual([
      { name: "cp1" },
      { name: "noncoding", sites: ["101-1600\\3", "7"], model: "jc69" },
    ]);
  });

  test("clearing an optional field leaves it out of the entry", async () => {
    render([{ name: "cp1", sites: ["1-9"], model: "infer" }]);
    const user = userEvent.setup({ document });

    await user.tripleClick(page().getByRole("textbox", { name: "Sites of Entries entry 1" }));
    await user.keyboard("{Backspace}");
    await user.selectOptions(page().getByRole("combobox", { name: "Model of Entries entry 1" }), "Default: infer");

    expect(settingValue()).toStrictEqual([{ name: "cp1" }]);
  });

  test("adding and removing entries changes the list", async () => {
    render([{ name: "cp1" }, { name: "cp2" }]);
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Add entry" }));
    const added = settingValue();
    await user.click(page().getByRole("button", { name: "Remove Entries entry 1" }));
    const removed = settingValue();

    const names = page()
      .getAllByRole<HTMLInputElement>("textbox", { name: /^Name/u })
      .map((field) => field.value);

    expect([added, removed, names]).toStrictEqual([
      [{ name: "cp1" }, { name: "cp2" }, {}],
      [{ name: "cp2" }, {}],
      ["cp2", ""],
    ]);
  });

  test("removing the last entry leaves an empty list", async () => {
    render([{ name: "cp1" }]);
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Remove Entries entry 1" }));

    expect([settingValue(), page().getByRole("grid", { name: "Entries" }).textContent]).toStrictEqual([
      [],
      "No entries",
    ]);
  });
});

function render(value: FormValue | undefined) {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(<Form value={value} />));
}

function settingValue(): unknown {
  return JSON.parse(page().getByLabelText("Entries setting").textContent ?? "null");
}

function page() {
  return within(document.body);
}

function Form({ value }: { value: FormValue | undefined }) {
  const form = useForm<FormConfig>({ defaultValues: { entries: value } });

  return (
    <FormProvider {...form}>
      <SettingControl command="ancestral" spec={ENTRIES} label="Entries" />
      <SettingValue />
    </FormProvider>
  );
}

function SettingValue() {
  const value = useWatch<FormConfig>({ name: "entries" });

  return <output aria-label="Entries setting">{JSON.stringify(value ?? null)}</output>;
}
