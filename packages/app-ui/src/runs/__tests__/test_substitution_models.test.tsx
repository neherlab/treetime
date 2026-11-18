import type { RunModels } from "@neherlab/app-contracts";
import { within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import { ModelsView } from "../SubstitutionModels";

const NUCLEOTIDES = ["A", "C", "G", "T"];

const UNIFORM = [
  [0, 1, 1, 1],
  [1, 0, 1, 1],
  [1, 1, 0, 1],
  [1, 1, 1, 0],
];

const MODELS: RunModels = {
  alignments: [{ alignment: 1, path: "/data/aln.fasta", sites: 900 }],
  models: [
    {
      kind: "sequence",
      partition: 1,
      name: "OPG001",
      alignment: 1,
      sites: 600,
      model: "infer",
      mu: 0.001,
      states: NUCLEOTIDES,
      pi: [0.25, 0.25, 0.25, 0.25],
      rates: UNIFORM,
      shared_params: ["exchange-rates"],
    },
    {
      kind: "sequence",
      partition: 2,
      name: "OPG002",
      alignment: 1,
      sites: 300,
      model: "infer",
      mu: 0.002,
      states: NUCLEOTIDES,
      pi: [0.25, 0.25, 0.25, 0.25],
      rates: UNIFORM,
      shared_params: ["exchange-rates"],
    },
    {
      kind: "trait",
      attribute: "country",
      model: "infer",
      mu: 0.5,
      states: ["UK", "US"],
      pi: [0.5, 0.5],
      rates: [
        [0, 0.5],
        [0.5, 0],
      ],
    },
  ],
  shared_params: ["exchange-rates"],
};

const HAPPY_DOM_LAYOUT_SIZES = ["clientWidth", "clientHeight"];

const mounted: Root[] = [];

beforeAll(() => {
  registerDom();
  renderAllVirtualRows();
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("substitution models", () => {
  test("selects the first model and shows its charts and shared parameters when the panel opens", () => {
    render();

    expect([
      selectedRows(),
      rateStates("Partition 1 (OPG001)"),
      frequencyChartShown("Partition 1 (OPG001)"),
      textShown("Partition 1 (OPG001)", "Shared parameters: exchange rates, shared by all 2 partitions"),
    ]).toStrictEqual([["Partition 1 (OPG001)"], ["A", "C", "G", "T"], true, true]);
  });

  test("shows the charts of the row that the user selects", async () => {
    render();
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("row", { name: /Trait country/u }));

    expect([selectedRows(), rateStates("Trait country")]).toStrictEqual([["Trait country"], ["UK", "US"]]);
  });
});

function render() {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(<ModelsView models={MODELS} />));
}

function selectedRows(): string[] {
  return page()
    .getAllByRole("row", { selected: true })
    .map((row) => within(row).getByRole("rowheader").textContent);
}

function frequencyChartShown(model: string): boolean {
  return (
    within(page().getByRole("region", { name: model })).queryByRole("figure", { name: "Equilibrium frequencies" }) !==
    null
  );
}

function textShown(model: string, text: string): boolean {
  return within(page().getByRole("region", { name: model })).queryByText(text) !== null;
}

function rateStates(model: string): string[] {
  const rates = within(page().getByRole("region", { name: model })).getByRole("figure", { name: "Substitution rates" });

  return within(rates)
    .getAllByRole("rowheader")
    .map((header) => header.textContent);
}

function page() {
  return within(document.body);
}

function renderAllVirtualRows() {
  for (const property of HAPPY_DOM_LAYOUT_SIZES) {
    Reflect.deleteProperty(HTMLElement.prototype, property);
  }
}
