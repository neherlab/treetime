import { auspiceEntries, transformAuspiceSource } from "../auspice-vite";

const AUSPICE_ROOT = "/repo/node_modules/.bun/auspice@3.1.0+ee1ed944166f6889/node_modules/auspice";

const TIMERS_SOURCE = `
export const timerStart = (name) => { console.log("start", name); };
export const timerEnd = (name) => { console.warn("end", name); };
`;

const RUNTIME_CONSOLE = Object.getOwnPropertyDescriptor(globalThis, "console");

describe("auspice source transform", () => {
  afterEach(() => {
    if (RUNTIME_CONSOLE !== undefined) {
      Object.defineProperty(globalThis, "console", RUNTIME_CONSOLE);
    }
  });

  test("the timing module becomes timers that print nothing", async () => {
    const printed: unknown[][] = [];

    const record = (...args: unknown[]) => {
      printed.push(args);
    };

    Object.defineProperty(globalThis, "console", {
      value: { log: record, warn: record, error: record },
      configurable: true,
    });

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

describe("auspice module entries", () => {
  test("each declaration names the Auspice module it declares", () => {
    expect(
      auspiceEntries([
        "reducers/controls.d.ts",
        "components/tree/index.d.ts",
        "components\\download\\helperFunctions.d.ts",
        "actions/types.d.ts",
      ]),
    ).toStrictEqual([
      "auspice/src/actions/types",
      "auspice/src/components/download/helperFunctions",
      "auspice/src/components/tree",
      "auspice/src/reducers/controls",
    ]);
  });

  test("folders and files other than declarations name no module", () => {
    expect(auspiceEntries(["components", "components/tree", "README.md", "reducers/tree/index.d.ts"])).toStrictEqual([
      "auspice/src/reducers/tree",
    ]);
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
