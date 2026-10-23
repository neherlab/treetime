import { within } from "@testing-library/dom";
import { userEvent } from "@testing-library/user-event";
import type { ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import { LocaleProvider, type Locale } from "../../i18n/locale";
import { ThemeProvider } from "../../ThemeProvider";
import { SearchTrigger } from "../SearchTrigger";
import { ThemeCycleButton, type ThemeChoice } from "../ThemeCycleButton";

const mounted: Root[] = [];

beforeAll(registerDom);

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
  localStorage.clear();
});

describe("theme cycle button", () => {
  test.each<[Locale, string[]]>([
    ["en", ["Theme follows the system", "Light theme", "Dark theme", "Theme follows the system"]],
    ["de", ["Design folgt dem System", "Helles Design", "Dunkles Design", "Design folgt dem System"]],
  ])("cycles system, light, and dark with their names in %s", async (locale, expected) => {
    render(
      <ThemeProvider>
        <ThemeCycleButton />
      </ThemeProvider>,
      locale,
    );
    const user = userEvent.setup({ document });
    const names = [buttonName()];

    for (const _ of expected.slice(1)) {
      await user.click(page().getByRole("button"));
      names.push(buttonName());
    }

    expect(names).toStrictEqual(expected);
  });

  test("reports each chosen theme to the surface", async () => {
    const chosen = new CallLog<ThemeChoice>();
    render(
      <ThemeProvider>
        <ThemeCycleButton onChosen={chosen.record} />
      </ThemeProvider>,
      "en",
    );
    const user = userEvent.setup({ document });

    await user.click(page().getByRole("button", { name: "Theme follows the system" }));
    await user.click(page().getByRole("button", { name: "Light theme" }));

    expect(chosen.calls).toStrictEqual(["light", "dark"]);
  });
});

describe("search trigger", () => {
  test("shows its hotkey and opens the search on click", async () => {
    const clicks = new CallLog<undefined>();
    render(<SearchTrigger label="Search the docs" hotkey="Mod+K" onClick={clicks.record} />, "en");
    const user = userEvent.setup({ document });
    const button = page().getByRole("button", { name: "Search the docs" });

    await user.click(button);

    expect([within(button).getByText(/K$/u).tagName, clicks.calls.length]).toStrictEqual(["KBD", 1]);
  });
});

function render(node: ReactNode, locale: Locale) {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);
  inAct(() => root.render(<LocaleProvider locale={locale}>{node}</LocaleProvider>));
}

class CallLog<T> {
  readonly calls: Array<T | undefined> = [];

  readonly record = (value?: T): void => {
    this.calls.push(value);
  };
}

function buttonName(): string | null {
  return page().getByRole("button").getAttribute("aria-label");
}

function page() {
  return within(document.body);
}
