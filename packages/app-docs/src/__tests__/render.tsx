import { inAct } from "@neherlab/app-ui/testing/dom";
import { within } from "@testing-library/dom";
import type { ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";

const mounted: Root[] = [];

export function setupDom(): void {
  afterEach(() => {
    for (const root of mounted.splice(0)) {
      inAct(() => root.unmount());
    }

    document.body.replaceChildren();
    localStorage.clear();
  });
}

export function render(node: ReactNode): void {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(node));
}

export function page() {
  return within(document.body);
}
