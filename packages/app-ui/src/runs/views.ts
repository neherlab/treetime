import type { AppCommand, RunStatus } from "@neherlab/app-contracts";

export type RunViewKey = "overview" | "clock" | "population" | "model" | "log" | "settings";

export interface RunViewTab {
  key: RunViewKey;
  label: string;
  to: RunViewRoute;
}

export type RunViewRoute = (typeof VIEW_TABS)[RunViewKey]["to"];

const VIEW_TABS = {
  overview: { label: "Overview", to: "/runs/$id/results" },
  clock: { label: "Clock", to: "/runs/$id/clock" },
  population: { label: "Population", to: "/runs/$id/population" },
  model: { label: "Model", to: "/runs/$id/model" },
  log: { label: "Log", to: "/runs/$id/log" },
  settings: { label: "Settings", to: "/runs/$id/settings" },
} as const satisfies Record<RunViewKey, { label: string; to: string }>;

const RESULT_VIEWS: Readonly<Record<AppCommand, readonly RunViewKey[]>> = {
  timetree: ["clock", "population", "model"],
  clock: [],
  ancestral: ["model"],
  optimize: ["model"],
  prune: [],
};

export function runViews(command: AppCommand, status: RunStatus): RunViewTab[] {
  const results = status === "ok" ? RESULT_VIEWS[command] : [];

  return (["overview", ...results, "log", "settings"] as const).map((key) => {
    const { label, to } = VIEW_TABS[key];

    return { key, label, to };
  });
}

export function hasRunView(command: AppCommand, status: RunStatus, key: RunViewKey): boolean {
  return runViews(command, status).some((view) => view.key === key);
}
