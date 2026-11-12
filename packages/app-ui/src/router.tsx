import { createRootRoute, createRoute, createRouter, Navigate } from "@tanstack/react-router";

import { DraftMain } from "./analysis/DraftMain";
import { WorkspaceRoute } from "./analysis/Workspace";
import { WORKSPACE_MAIN_ID } from "./analysis/WorkspaceLayout";
import { NotFoundPage } from "./NotFoundPage";
import { MAIN_SCROLL_ID, RootLayout } from "./RootLayout";
import { ComparePage } from "./runs/ComparePage";
import { RunLogTab, RunPage, RunResultsTab, RunSettingsTab } from "./runs/RunPage";

const rootRoute = createRootRoute({ component: RootLayout });

const indexRoute = createRoute({ getParentRoute: () => rootRoute, path: "/", component: IndexRedirect });

const workspaceRoute = createRoute({ getParentRoute: () => rootRoute, id: "workspace", component: WorkspaceRoute });

const newRoute = createRoute({ getParentRoute: () => workspaceRoute, path: "/new", component: DraftMain });

const runRoute = createRoute({ getParentRoute: () => workspaceRoute, path: "/runs/$id", component: RunLayout });

const runIndexRoute = createRoute({ getParentRoute: () => runRoute, path: "/", component: RunRedirect });

const runResultsRoute = createRoute({ getParentRoute: () => runRoute, path: "/results", component: RunResultsTab });

const runSettingsRoute = createRoute({ getParentRoute: () => runRoute, path: "/settings", component: RunSettingsTab });

const runLogRoute = createRoute({ getParentRoute: () => runRoute, path: "/log", component: RunLogTab });

const compareRoute = createRoute({ getParentRoute: () => rootRoute, path: "/compare/$a/$b", component: CompareRuns });

const routeTree = rootRoute.addChildren([
  indexRoute,
  workspaceRoute.addChildren([
    newRoute,
    runRoute.addChildren([runIndexRoute, runResultsRoute, runSettingsRoute, runLogRoute]),
  ]),
  compareRoute,
]);

export const router = createRouter({
  routeTree,
  defaultPreload: "intent",
  defaultNotFoundComponent: NotFoundPage,
  scrollRestoration: true,
  scrollToTopSelectors: [`#${MAIN_SCROLL_ID}`, `#${WORKSPACE_MAIN_ID}`],
});

declare module "@tanstack/react-router" {
  interface Register {
    router: typeof router;
  }
}

function IndexRedirect() {
  return <Navigate to="/new" replace />;
}

function RunLayout() {
  const { id } = runRoute.useParams();

  return <RunPage id={id} />;
}

function RunRedirect() {
  const { id } = runRoute.useParams();

  return <Navigate to="/runs/$id/results" params={{ id }} replace />;
}

function CompareRuns() {
  const { a, b } = compareRoute.useParams();

  return <ComparePage first={a} second={b} />;
}
