import { type ApiClient, createApiClient } from "@neherlab/app-contracts/client";

export function createWebApiClient(): ApiClient {
  return createApiClient({ baseUrl: globalThis.location.origin, fetch: globalThis.fetch.bind(globalThis) });
}
