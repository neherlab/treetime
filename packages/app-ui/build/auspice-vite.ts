import { transformWithOxc, type Plugin } from "vite";

const AUSPICE_SOURCE =
  /[\\/]node_modules[\\/](?:\.bun[\\/][^\\/]+[\\/]node_modules[\\/])?auspice[\\/]src[\\/][^?]+\.js(?:\?.*)?$/u;

const AUSPICE_ENTRIES = [
  "auspice/src/actions/colors",
  "auspice/src/actions/recomputeReduxState",
  "auspice/src/actions/tree",
  "auspice/src/actions/types",
  "auspice/src/components/controls/choose-branch-labelling",
  "auspice/src/components/controls/choose-layout",
  "auspice/src/components/controls/choose-metric",
  "auspice/src/components/controls/choose-tip-label",
  "auspice/src/components/controls/color-by",
  "auspice/src/components/controls/controlHeader",
  "auspice/src/components/controls/filter",
  "auspice/src/components/controls/miscInfoText",
  "auspice/src/components/controls/styles",
  "auspice/src/components/controls/toggle-focus",
  "auspice/src/components/download/downloadButtons",
  "auspice/src/components/download/downloadModal",
  "auspice/src/components/info/filtersSummary",
  "auspice/src/components/tree",
  "auspice/src/middleware/performanceFlags",
  "auspice/src/middleware/scatterplot",
  "auspice/src/reducers/browserDimensions",
  "auspice/src/reducers/controls",
  "auspice/src/reducers/entropy",
  "auspice/src/reducers/frequencies",
  "auspice/src/reducers/measurements",
  "auspice/src/reducers/metadata",
  "auspice/src/reducers/narrative",
  "auspice/src/reducers/notifications",
  "auspice/src/reducers/tree",
  "auspice/src/reducers/tree/treeToo",
  "auspice/src/util/computeResponsive",
];

const APP_UI_PACKAGE = "@neherlab/app-ui";

const SHARED_PACKAGES = ["react", "react-dom"];

const AUSPICE_TIMERS = /[\\/]auspice[\\/]src[\\/]util[\\/]perf\.js(?:\?.*)?$/u;

const AUSPICE_TIMERS_STUB = "export const timerStart = () => {};\nexport const timerEnd = () => {};\n";

const AUSPICE_TRANSFORM = {
  filter: { id: AUSPICE_SOURCE },
  handler: transformAuspiceSource,
} satisfies Plugin["transform"];

export async function transformAuspiceSource(code: string, id: string) {
  if (AUSPICE_TIMERS.test(id)) {
    return { code: AUSPICE_TIMERS_STUB, map: null };
  }

  const result = await transformWithOxc(code, id, {
    lang: "jsx",
    jsx: { runtime: "automatic" },
    decorator: { legacy: true },
    sourcemap: true,
  });

  return { code: result.code, map: result.map ?? null };
}

export function auspice(): Plugin {
  return {
    name: "treetime-auspice",
    enforce: "pre",
    config() {
      return {
        define: {
          "process.env.EXTENSION_DATA": JSON.stringify(""),
          "process.env.SKIP_REDUX_CHECKS": JSON.stringify(""),
          "process.env.ENABLE_SERVICE_WORKER": "false",
        },
        resolve: { dedupe: SHARED_PACKAGES },
        optimizeDeps: {
          extensions: [".tsx"],
          include: AUSPICE_ENTRIES.map((entry) => `${APP_UI_PACKAGE} > ${entry}`),
          rolldownOptions: { plugins: [{ name: "treetime-auspice-source", transform: AUSPICE_TRANSFORM }] },
        },
      };
    },
    transform: AUSPICE_TRANSFORM,
  };
}
