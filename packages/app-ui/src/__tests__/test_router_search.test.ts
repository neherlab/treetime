import { createMemoryHistory, createRootRoute, createRoute, createRouter } from "@tanstack/react-router";

import { rootSearch } from "../routerSearch";

describe("router search parameters", () => {
  test("a link keeps the jobs parameter of the current URL", async () => {
    const router = await testRouter("/?jobs=2");

    expect(router.buildLocation({ to: "/runs/$id", params: { id: "r1" } }).href).toBe("/runs/r1?jobs=2");
  });

  test("a jobs value that is not a number stays as written", async () => {
    const router = await testRouter("/runs/r1?jobs=many");

    expect(router.buildLocation({ to: "/new" }).href).toBe("/new?jobs=many");
  });

  test("a link keeps no other parameter of the current URL", async () => {
    const router = await testRouter("/runs/r1?jobs=2&view=tree");

    expect(router.buildLocation({ to: "/new" }).href).toBe("/new?jobs=2");
  });

  test("a URL without jobs gains none", async () => {
    const router = await testRouter("/");

    expect(router.buildLocation({ to: "/new" }).href).toBe("/new");
  });
});

async function testRouter(url: string) {
  const rootRoute = createRootRoute(rootSearch);

  const routeTree = rootRoute.addChildren([
    createRoute({ getParentRoute: () => rootRoute, path: "/" }),
    createRoute({ getParentRoute: () => rootRoute, path: "/new" }),
    createRoute({ getParentRoute: () => rootRoute, path: "/runs/$id" }),
  ]);

  const router = createRouter({ routeTree, history: createMemoryHistory({ initialEntries: [url] }) });

  await router.load();

  return router;
}
