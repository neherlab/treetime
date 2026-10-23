import "../../__tests__/register-dom";
import { inAct } from "@neherlab/app-ui/testing/dom";
import { userEvent } from "@testing-library/user-event";

import { page, render, setupDom } from "../../__tests__/render";
import type { NavNode } from "../../lib/nav";
import { CalloutContainer, CalloutDescription, CalloutTitle } from "../Callout";
import { Disclosure } from "../Disclosure";
import { FallbackNotice } from "../FallbackNotice";
import { NavTree } from "../NavTree";
import { SyncedTabs } from "../SyncedTabs";
import { TableOfContents } from "../TableOfContents";

setupDom();

describe("callouts", () => {
  test.each([
    ["info", "default"],
    ["success", "success"],
    ["warning", "warning"],
    ["error", "destructive"],
  ])("a %s callout is an app alert with the role of a note", (type, variant) => {
    render(
      <CalloutContainer type={type}>
        <CalloutTitle>Title</CalloutTitle>
        <CalloutDescription>Text</CalloutDescription>
      </CalloutContainer>,
    );
    const note = page().getByRole("note");

    expect([note.dataset["slot"], note.dataset["variant"], note.textContent]).toStrictEqual([
      "alert",
      variant,
      "TitleText",
    ]);
  });
});

describe("navigation tree", () => {
  test("marks the current page and opens the folder that holds it", () => {
    render(<NavTree nodes={NODES} currentUrl="/dev/reference/cli" locale="en" />);
    const current = page().getByRole("link", { current: "page" });

    const folders = [...document.querySelectorAll("details")].map((folder) => [
      folder.querySelector("summary")?.textContent,
      folder.open,
    ]);

    expect([current.textContent, folders]).toStrictEqual([
      "CLI",
      [
        ["Showcase", false],
        ["Reference", true],
      ],
    ]);
  });

  test("names the navigation in the page locale", () => {
    render(<NavTree nodes={NODES} currentUrl="/dev" locale="de" />);

    expect(page().getByRole("navigation").getAttribute("aria-label")).toBe("Doku-Navigation");
  });
});

describe("table of contents", () => {
  test("marks the heading that the anchor observer reports as visible", () => {
    const observers = installIntersectionObserver();
    render(
      <>
        <h2 id="first">First</h2>
        <h2 id="second">Second</h2>
        <TableOfContents headings={HEADINGS} locale="en" />
      </>,
    );

    inAct(() => observers.report("second"));

    expect(page().getByRole("link", { current: "location" }).textContent).toBe("Second");
  });
});

describe("fallback notice", () => {
  test("tells the reader in German that the page shows English text", () => {
    render(<FallbackNotice locale="de" />);

    expect(page().getByRole("note").textContent).toContain("Diese Seite ist noch nicht übersetzt");
  });
});

describe("synced tabs", () => {
  test("a choice in one group switches every group with the same key", async () => {
    render(
      <>
        <SyncedTabs syncKey="surface" labels={TAB_LABELS} tab0="first cli" tab1="first web" />
        <SyncedTabs syncKey="surface" labels={TAB_LABELS} tab0="second cli" tab1="second web" />
      </>,
    );
    const user = userEvent.setup({ document });
    const [firstWebTab] = page().getAllByRole("tab", { name: "Web app" });

    await user.click(firstWebTab ?? document.body);

    expect(
      page()
        .getAllByRole("tabpanel")
        .map((panel) => panel.textContent),
    ).toStrictEqual(["first web", "second web"]);
  });
});

describe("disclosure", () => {
  test("shows its content when opened", async () => {
    render(<Disclosure title="Question">Answer</Disclosure>);
    const user = userEvent.setup({ document });
    const answer = page().getByText("Answer");
    const hiddenBefore = answer.closest("[hidden]") !== null;

    await user.click(page().getByRole("button", { name: "Question" }));

    expect([hiddenBefore, answer.closest("[hidden]") !== null]).toStrictEqual([true, false]);
  });
});

const HEADINGS = [
  { depth: 2, slug: "first", text: "First" },
  { depth: 2, slug: "second", text: "Second" },
];

const TAB_LABELS = ["CLI", "Web app"];

const NODES: NavNode[] = [
  { kind: "page", name: "Home", url: "/dev" },
  {
    kind: "folder",
    name: "Showcase",
    url: undefined,
    children: [{ kind: "page", name: "Math", url: "/dev/showcase/math" }],
  },
  {
    kind: "folder",
    name: "Reference",
    url: undefined,
    children: [{ kind: "page", name: "CLI", url: "/dev/reference/cli" }],
  },
];

function installIntersectionObserver(): ObserverReports {
  const callbacks: ObserverCallback[] = [];

  class FakeIntersectionObserver {
    constructor(callback: ObserverCallback) {
      callbacks.push(callback);
    }

    observe(): void {}

    disconnect(): void {}
  }

  Object.assign(globalThis, { IntersectionObserver: FakeIntersectionObserver });

  return {
    report: (id) => {
      for (const callback of callbacks) {
        callback([{ target: document.getElementById(id), isIntersecting: true }]);
      }
    },
  };
}

interface ObserverReports {
  report: (id: string) => void;
}

type ObserverCallback = (entries: Array<{ target: Element | null; isIntersecting: boolean }>) => void;
