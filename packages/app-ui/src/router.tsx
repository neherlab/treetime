import { createRootRoute, createRoute, createRouter, Navigate } from "@tanstack/react-router";

import { NewAnalysisPage } from "./analysis/NewAnalysisPage";
import { RootLayout } from "./RootLayout";
import { ComparePage } from "./runs/ComparePage";
import { RunPage } from "./runs/RunPage";

const rootRoute = createRootRoute({ component: RootLayout });

const indexRoute = createRoute({ getParentRoute: () => rootRoute, path: "/", component: IndexRedirect });

const newRoute = createRoute({ getParentRoute: () => rootRoute, path: "/new", component: NewAnalysisPage });

const runRoute = createRoute({ getParentRoute: () => rootRoute, path: "/runs/$id", component: RunRedirect });

const runResultsRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "/runs/$id/results",
  component: RunResults,
});

const runSettingsRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "/runs/$id/settings",
  component: RunSettings,
});

const runLogRoute = createRoute({ getParentRoute: () => rootRoute, path: "/runs/$id/log", component: RunLog });

const compareRoute = createRoute({ getParentRoute: () => rootRoute, path: "/compare/$a/$b", component: CompareRuns });

const routeTree = rootRoute.addChildren([
  indexRoute,
  newRoute,
  runRoute,
  runResultsRoute,
  runSettingsRoute,
  runLogRoute,
  compareRoute,
]);

export const router = createRouter({ routeTree, defaultPreload: "intent" });

declare module "@tanstack/react-router" {
  interface Register {
    router: typeof router;
  }
}

function IndexRedirect() {
  return <Navigate to="/new" replace />;
}

function RunRedirect() {
  const { id } = runRoute.useParams();

  return <Navigate to="/runs/$id/results" params={{ id }} replace />;
}

function RunResults() {
  const { id } = runResultsRoute.useParams();

  return <RunPage id={id} tab="results" />;
}

function RunSettings() {
  const { id } = runSettingsRoute.useParams();

  return <RunPage id={id} tab="settings" />;
}

function RunLog() {
  const { id } = runLogRoute.useParams();

  return <RunPage id={id} tab="log" />;
}

function CompareRuns() {
  const { a, b } = compareRoute.useParams();

  return <ComparePage first={a} second={b} />;
}
