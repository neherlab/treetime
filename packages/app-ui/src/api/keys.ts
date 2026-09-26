import type { StalePath, StaleScope } from "@neherlab/app-contracts";
import type { ApiClient, ApiRequestOptions } from "@neherlab/app-contracts/client";

export interface ApiCallContext {
  client: ApiClient;
  throwOnError: true;
  signal?: AbortSignal;
}

export type ApiRequest = (context: ApiCallContext) => Promise<object>;

export type ApiKey = readonly unknown[];

type CapturedRequest = Pick<ApiRequestOptions, "url" | "path" | "query" | "body">;

type RequestQuery = CapturedRequest["query"];

const SCOPE_COVERS_LENGTH: Record<StaleScope, (keyPathLength: number, stalePathLength: number) => boolean> = {
  exact: (keyPathLength, stalePathLength) => keyPathLength === stalePathLength,
  subtree: (keyPathLength, stalePathLength) => keyPathLength >= stalePathLength,
};

export function staleCoversKey(stale: StalePath, key: readonly unknown[]): boolean {
  const stalePath = pathKey(stale.path);

  return (
    SCOPE_COVERS_LENGTH[stale.scope](pathLengthOfKey(key), stalePath.length) &&
    stalePath.every((segment, index) => key[index] === segment)
  );
}

export function requestKey(client: ApiClient, request: ApiRequest): ApiKey {
  const captured = captureRequest(client, request);
  const path = client.buildUrl({ url: captured.url, path: captured.path ?? {}, baseUrl: "" });
  const query = sortedQuery(captured.query);

  return [
    ...pathKey(path),
    ...(query === undefined ? [] : [query]),
    ...(captured.body === undefined ? [] : [captured.body]),
  ];
}

function pathKey(path: string): string[] {
  return path
    .split("/")
    .filter((segment) => segment !== "")
    .map((segment) => decodeURIComponent(segment));
}

function pathLengthOfKey(key: readonly unknown[]): number {
  const end = key.findIndex((part) => !isPathSegment(part));

  return end === -1 ? key.length : end;
}

function isPathSegment(part: unknown): part is string {
  return typeof part === "string";
}

function captureRequest(client: ApiClient, request: ApiRequest): CapturedRequest {
  const requests: CapturedRequest[] = [];

  const record = (options: CapturedRequest): Promise<never> => {
    requests.push(options);

    return Promise.withResolvers<never>().promise;
  };

  const capturing: ApiClient = {
    ...client,
    connect: record,
    delete: record,
    get: record,
    head: record,
    options: record,
    patch: record,
    post: record,
    put: record,
    request: record,
    trace: record,
    sse: {
      connect: record,
      delete: record,
      get: record,
      head: record,
      options: record,
      patch: record,
      post: record,
      put: record,
      trace: record,
    },
  };

  void request({ client: capturing, throwOnError: true });

  const [captured, ...rest] = requests;

  if (captured === undefined || rest.length > 0) {
    throw new TypeError(`the request sent ${requests.length} requests to its client instead of one`);
  }

  return captured;
}

function sortedQuery(query: RequestQuery): RequestQuery {
  const entries = Object.entries(query ?? {})
    .filter(([, value]) => value !== undefined)
    .toSorted(([a], [b]) => a.localeCompare(b));

  return entries.length === 0 ? undefined : Object.fromEntries(entries);
}
