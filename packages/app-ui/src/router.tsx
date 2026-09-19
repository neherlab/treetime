import { createRootRoute, createRoute, createRouter, Navigate } from "@tanstack/react-router";

import { Workspace } from "./components/Workspace";
import { RootLayout } from "./RootLayout";
import { DEFAULT_COMMAND } from "./types";

const rootRoute = createRootRoute({ component: RootLayout });

function IndexRedirect() {
  return <Navigate to="/$command" params={{ command: DEFAULT_COMMAND }} replace />;
}

const indexRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "/",
  component: IndexRedirect,
});

const commandRoute = createRoute({
  getParentRoute: () => rootRoute,
  path: "$command",
  component: Workspace,
});

const routeTree = rootRoute.addChildren([indexRoute, commandRoute]);

export const router = createRouter({
  routeTree,
  defaultPreload: "intent",
});

declare module "@tanstack/react-router" {
  interface Register {
    router: typeof router;
  }
}
