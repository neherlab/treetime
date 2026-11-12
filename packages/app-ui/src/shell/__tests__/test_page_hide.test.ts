import { onPageHide } from "../pageHide";

describe("work that a page hide commits", () => {
  test("runs once when the page hides", () => {
    const page = new EventTarget();
    const calls: string[] = [];
    onPageHide(page, () => {
      calls.push("commit");
    });

    page.dispatchEvent(new Event("pagehide"));
    page.dispatchEvent(new Event("pagehide"));

    expect(calls).toStrictEqual(["commit"]);
  });

  test("does not run after it was stopped, for example by an undo", () => {
    const page = new EventTarget();
    const calls: string[] = [];

    const stop = onPageHide(page, () => {
      calls.push("commit");
    });

    stop();
    page.dispatchEvent(new Event("pagehide"));

    expect(calls).toStrictEqual([]);
  });
});
