import type { SettingSpec } from "@neherlab/app-contracts";
import { within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import { createRoot, type Root } from "react-dom/client";
import { FormProvider, useForm, useWatch } from "react-hook-form";

import { inAct, registerDom } from "../../__tests__/dom";
import type { FormConfig, FormValue } from "../formValues";
import { SettingControl } from "../SettingControl";

const COALESCENT: SettingSpec = {
  key: "coalescent",
  label: "Coalescent prior",
  flag: "--coalescent",
  group: "Coalescent prior",
  role: "setting",
  kind: "enum-or-text",
  nullable: true,
  options: [
    { value: "opt", help: "" },
    { value: "skyline", help: "" },
  ],
  item_kind: "string",
  fields: [],
  examples: [0.01],
  value_names: ["TC|opt|skyline"],
  conflicts: [],
  help: "Coalescent prior",
  more: "",
};

const MODEL: SettingSpec = {
  key: "model",
  label: "Substitution model",
  flag: "--model",
  group: "Substitution model",
  role: "setting",
  kind: "enum",
  nullable: true,
  options: [
    { value: "infer", help: "" },
    { value: "jc69", help: "" },
  ],
  item_kind: "string",
  fields: [],
  documented_default: "infer",
  examples: [],
  value_names: ["MODEL"],
  conflicts: [],
  help: "Substitution model to use",
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

describe("coalescent prior control", () => {
  test.each([
    ["an unset prior", undefined, "Not set", null],
    ["a time scale in substitutions per site", 0.01, "Value", "0.01"],
    ["a time scale with a unit", "10y", "Value", "10y"],
    ["a fitted constant time scale", "opt", "opt", null],
    ["a skyline", "skyline", "skyline", null],
  ])("shows %s", (_, value, option, text) => {
    render(value);

    expect([selectedOption(), valueField()?.value ?? null]).toStrictEqual([option, text]);
  });

  test("picking an option, a value and no prior writes each to the setting", async () => {
    render(undefined);
    const user = userEvent.setup({ document });

    await user.selectOptions(page().getByRole("combobox", { name: "Coalescent prior" }), "skyline");
    const skyline = settingText();
    await user.selectOptions(page().getByRole("combobox", { name: "Coalescent prior" }), "Value");
    const example = [settingText(), valueField()?.value];
    const field = page().getByRole("textbox", { name: "Coalescent prior value" });
    await user.clear(field);
    await user.type(field, "26w");
    const typed = settingText();
    await user.selectOptions(page().getByRole("combobox", { name: "Coalescent prior" }), "Not set");

    expect([skyline, example, typed, settingText()]).toStrictEqual(['"skyline"', ["0.01", "0.01"], '"26w"', "null"]);
  });
});

describe("setting with a documented default", () => {
  test("an unset setting names its documented default and a pick of that value sets it", async () => {
    render(undefined, MODEL);
    const user = userEvent.setup({ document });
    const select = page().getByRole<HTMLSelectElement>("combobox", { name: "Substitution model" });
    const unset = [select.selectedOptions[0]?.textContent, settingText()];

    await user.selectOptions(select, "infer");

    expect([unset, settingText()]).toStrictEqual([["Default: infer", "null"], '"infer"']);
  });
});

function render(value: FormValue | undefined, spec: SettingSpec = COALESCENT) {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(<Form spec={spec} value={value} />));
}

function selectedOption(): string | undefined {
  const select = page().getByRole<HTMLSelectElement>("combobox", { name: "Coalescent prior" });

  return select.selectedOptions[0]?.textContent ?? undefined;
}

function valueField(): HTMLInputElement | null {
  return page().queryByRole<HTMLInputElement>("textbox", { name: "Coalescent prior value" });
}

function settingText(): string | null {
  return page().getByLabelText("Setting value").textContent;
}

function page() {
  return within(document.body);
}

function Form({ spec, value }: { spec: SettingSpec; value: FormValue | undefined }) {
  const form = useForm<FormConfig>({ defaultValues: { [spec.key]: value } });

  return (
    <FormProvider {...form}>
      <SettingControl command="timetree" spec={spec} label={spec.label} />
      <SettingValue name={spec.key} />
    </FormProvider>
  );
}

function SettingValue({ name }: { name: string }) {
  const value = useWatch<FormConfig>({ name });

  return <output aria-label="Setting value">{JSON.stringify(value ?? null)}</output>;
}
