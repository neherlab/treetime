import { within } from "@testing-library/dom";
import type { ReactElement } from "react";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import type { ViewPanel } from "../layout";
import { PanelBoundary } from "../PanelBoundary";

const mounted: Root[] = [];

const failures: ViewPanel[] = [];

beforeAll(() => {
  registerDom();
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  failures.splice(0);
  document.body.replaceChildren();
});

describe("panel error boundary", () => {
  test("a panel that fails shows its error in its place and leaves the other panels", () => {
    quietly(() =>
      render(
        <>
          <PanelBoundary panel="tree" onFailure={recordFailure}>
            <p>Tree drawing</p>
          </PanelBoundary>
          <PanelBoundary panel="entropy" onFailure={recordFailure}>
            <FailingPanel message="Cannot read properties of undefined (reading 'bars')" />
          </PanelBoundary>
        </>,
      ),
    );

    const page = within(document.body);

    expect([page.queryByText("Tree drawing") !== null, page.getByRole("alert").textContent, failures]).toStrictEqual([
      true,
      "Auspice cannot draw the entropy chartCannot read properties of undefined (reading 'bars')",
      ["entropy"],
    ]);
  });
});

function recordFailure(panel: ViewPanel): void {
  failures.push(panel);
}

function FailingPanel({ message }: { message: string }): never {
  throw new Error(message);
}

function render(element: ReactElement): void {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);

  inAct(() => root.render(element));
}

function quietly(run: () => void): void {
  const original = console.error;
  console.error = () => undefined;

  try {
    run();
  } finally {
    console.error = original;
  }
}
