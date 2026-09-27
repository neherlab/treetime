import { afterEach, describe, expect, test, vi } from "vitest";

import { transformAuspiceSource } from "../auspice-vite";

const AUSPICE_ROOT = "/repo/node_modules/.bun/auspice@3.0.0+ee1ed944166f6889/node_modules/auspice";

const TIMERS_SOURCE = `
export const timerStart = (name) => { console.log("start", name); };
export const timerEnd = (name) => { console.warn("end", name); };
`;

describe("auspice source transform", () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  test("the timing module becomes timers that print nothing", async () => {
    const printed: unknown[][] = [];

    const record = (...args: unknown[]) => {
      printed.push(args);
    };

    vi.stubGlobal("console", { log: record, warn: record, error: record });

    const { code } = await transformAuspiceSource(TIMERS_SOURCE, `${AUSPICE_ROOT}/src/util/perf.js`);
    const timers = await importModule(code);

    timers.timerStart("modifySVG");
    timers.timerStart("modifySVG");
    timers.timerEnd("modifySVG");
    timers.timerEnd("mapToScreen");

    expect(printed).toStrictEqual([]);
  });

  test("the timing module requested with a query suffix is also replaced", async () => {
    const { code } = await transformAuspiceSource(TIMERS_SOURCE, `${AUSPICE_ROOT}/src/util/perf.js?v=1a2b3c`);

    expect(code).not.toContain("console");
  });

  test("a module named like the timing module elsewhere in auspice keeps its code", async () => {
    const { code } = await transformAuspiceSource(TIMERS_SOURCE, `${AUSPICE_ROOT}/src/components/perf.js`);

    expect(code).toContain('console.log("start", name)');
  });

  test("other auspice modules have their JSX compiled", async () => {
    const { code } = await transformAuspiceSource(
      "export const Label = () => <span>tree</span>;",
      `${AUSPICE_ROOT}/src/components/label.js`,
    );

    expect(code).toContain("react/jsx-runtime");
    expect(code).not.toContain("<span>");
  });
});

async function importModule(code: string): Promise<Timers> {
  const module: unknown = await import(`data:text/javascript,${encodeURIComponent(code)}`);

  if (!isTimers(module)) {
    throw new Error("the transformed timing module does not export timerStart and timerEnd");
  }

  return module;
}

function isTimers(module: unknown): module is Timers {
  return (
    typeof module === "object" &&
    module !== null &&
    "timerStart" in module &&
    typeof module.timerStart === "function" &&
    "timerEnd" in module &&
    typeof module.timerEnd === "function"
  );
}

interface Timers {
  timerStart(name: string): void;
  timerEnd(name: string): void;
}
