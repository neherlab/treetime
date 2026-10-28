import { within } from "@testing-library/dom";
import type { ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import { createAuspiceStore, type AuspiceStore } from "../store";

const mounted: Root[] = [];

beforeAll(() => {
  registerDom();
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("auspice messages", () => {
  test("a message shows as a toast and leaves Auspice's message list", async () => {
    const store = createAuspiceStore();
    await render(store);

    inAct(() =>
      store.dispatch({
        type: "ADD_NOTIFICATION",
        id: 1,
        notificationType: "error",
        message: "Error saving tree!",
        details: "An empty tree was created.",
      }),
    );

    const page = within(document.body);

    expect([
      await page.findByText("Error saving tree!"),
      page.getByText("An empty tree was created."),
      store.getState().notifications.stack,
    ]).toStrictEqual([expect.any(HTMLElement), expect.any(HTMLElement), []]);
  });

  test("messages from before the viewer opened show once each", async () => {
    const store = createAuspiceStore();

    store.dispatch({ type: "ADD_NOTIFICATION", id: 1, notificationType: "info", message: "Tree written", details: "" });
    store.dispatch({
      type: "ADD_NOTIFICATION",
      id: 2,
      notificationType: "warning",
      message: "Map skipped",
      details: "",
    });
    await render(store);

    const page = within(document.body);

    expect([
      await page.findAllByText("Tree written"),
      await page.findAllByText("Map skipped"),
      store.getState().notifications.stack,
    ]).toStrictEqual([[expect.any(HTMLElement)], [expect.any(HTMLElement)], []]);
  });
});

async function render(store: AuspiceStore) {
  const { Toaster } = await import("@neherlab/design/ui/toast");
  const { useNotificationToasts } = await import("../notifications");

  function Messages(): ReactNode {
    useNotificationToasts(store);

    return null;
  }

  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() =>
    root.render(
      <Toaster>
        <Messages />
      </Toaster>,
    ),
  );
}
