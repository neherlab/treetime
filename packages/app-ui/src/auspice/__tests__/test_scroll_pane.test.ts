import { registerDom } from "../../__tests__/dom";
import { scrollAncestor } from "../scrollPane";

describe("scroll pane of the tree viewer", () => {
  beforeAll(() => {
    registerDom();
  });

  test("is the nearest ancestor that scrolls, not the outer page scroll", () => {
    const page = nested(document.body, "auto");
    const pane = nested(page, "auto");
    const viewer = nested(nested(pane, "visible"), "visible");

    expect(scrollAncestor(viewer)).toBe(pane);
  });

  test("is absent when no ancestor scrolls", () => {
    const viewer = nested(nested(document.body, "visible"), "visible");

    expect(scrollAncestor(viewer)).toBeUndefined();
  });
});

function nested(parent: HTMLElement, overflowY: string): HTMLElement {
  const element = document.createElement("div");
  element.style.overflowY = overflowY;
  parent.append(element);

  return element;
}
