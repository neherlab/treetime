import {
  events,
  resumableStream,
  runsEvents,
  SSE_MAX_RETRY_ATTEMPTS,
  type ApiClient,
  type ResumableStreamOptions,
} from "@neherlab/app-contracts/client";
import {
  experimental_streamedQuery,
  queryOptions,
  useQuery,
  useQueryClient,
  type QueryClient,
  type UseQueryResult,
} from "@tanstack/react-query";
import { useEffect, useState } from "react";

import { EMPTY_PROGRESS, foldRunEvents, type RunEvent, type RunProgress } from "../results/progress";
import { useApiContext } from "./context";
import { pathKey, requestKey } from "./keys";

type StreamTiming = Pick<ResumableStreamOptions<{ seq: number }>, "retryDelay" | "maxRetryDelay" | "sleep">;

export interface FollowAppEventsOptions extends StreamTiming {
  client: ApiClient;
  queryClient: QueryClient;
  signal: AbortSignal;
}

export interface RunEventStreamOptions extends StreamTiming {
  from?: number | undefined;
  signal?: AbortSignal | undefined;
}

export function useAppEvents(): Error | undefined {
  const { client } = useApiContext();
  const queryClient = useQueryClient();
  const [failure, setFailure] = useState<Error | undefined>(undefined);

  useEffect(() => {
    const controller = new AbortController();

    const follow = async () => {
      try {
        await followAppEvents({ client, queryClient, signal: controller.signal });
      } catch (error: unknown) {
        if (!controller.signal.aborted) {
          const failed = error instanceof Error ? error : new Error(String(error));
          console.error("[TreeTime] the app event stream failed", failed);
          setFailure(failed);
        }
      }
    };

    void follow();

    return () => controller.abort();
  }, [client, queryClient]);

  return failure;
}

export async function followAppEvents({
  client,
  queryClient,
  signal,
  ...timing
}: FollowAppEventsOptions): Promise<void> {
  const stream = resumableStream({ ...timing, open: (request) => events({ ...request, client }), from: 0, signal });

  for await (const event of stream) {
    void invalidateStale(queryClient, event.stale);
  }
}

export async function invalidateStale(queryClient: QueryClient, stale: readonly string[]): Promise<void> {
  await Promise.all(stale.map((path) => queryClient.invalidateQueries({ queryKey: pathKey(path) })));
}

export function useRunEvents(id: string): UseQueryResult<RunProgress> {
  const { client } = useApiContext();

  return useQuery(runEventsQueryOptions(client, id));
}

export function runEventsQueryOptions(client: ApiClient, id: string, timing: StreamTiming = {}) {
  return queryOptions<RunProgress>({
    queryKey: requestKey(client, (context) => runsEvents({ ...context, path: { id } })),
    queryFn: experimental_streamedQuery<RunEvent, RunProgress>({
      streamFn: ({ signal }) => runEventStream(client, id, { ...timing, signal }),
      reducer: (progress, event) => foldRunEvents(progress, [event]),
      initialValue: EMPTY_PROGRESS,
    }),
    staleTime: "static",
    retry: false,
  });
}

export async function* runEventStream(
  client: ApiClient,
  id: string,
  { from = 0, signal, ...timing }: RunEventStreamOptions = {},
): AsyncGenerator<RunEvent> {
  const stream = resumableStream({
    ...timing,
    open: (request) => runsEvents({ ...request, client, path: { id } }),
    from,
    signal,
    maxAttempts: SSE_MAX_RETRY_ATTEMPTS,
  });

  for await (const event of stream) {
    yield event;

    if (event.type === "terminal") {
      return;
    }
  }
}
