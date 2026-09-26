import { ZodError } from "zod";

import { type Client, type Config, createClient, createConfig, type RequestOptions } from "./generated/client";
import type { ErrorResponse } from "./generated/types.gen";
import { zErrorResponse } from "./generated/zod.gen";

export * from "./generated/sdk.gen";

export type { Client as ApiClient, RequestOptions as ApiRequestOptions } from "./generated/client";

export const SSE_MAX_RETRY_ATTEMPTS = 8;

export const SSE_RETRY_DELAY_MS = 1000;

export const SSE_MAX_RETRY_DELAY_MS = 30_000;

export interface ApiClientOptions {
  baseUrl: string;
  fetch: typeof globalThis.fetch;
}

export class ApiError extends Error {
  readonly status: number;
  readonly response: ErrorResponse;

  constructor(status: number, response: ErrorResponse) {
    super([response.message, ...response.causes].join(": "));
    this.name = "ApiError";
    this.status = status;
    this.response = response;
  }
}

export function createApiClient(options: ApiClientOptions): Client {
  const config: Config & Pick<RequestOptions, "onSseError" | "sseMaxRetryAttempts"> = {
    ...createConfig({ baseUrl: options.baseUrl, fetch: rejectErrorResponses(options.fetch) }),
    onSseError: (error) => {
      if (isFatalStreamError(error)) {
        throw error;
      }
    },
    sseMaxRetryAttempts: SSE_MAX_RETRY_ATTEMPTS,
  };

  return createClient(config);
}

export interface ResumableStreamRequest {
  query: { from?: number };
  signal: AbortSignal;
  sseMaxRetryAttempts: number;
  onSseError: NonNullable<RequestOptions["onSseError"]>;
}

export interface ResumableStreamOptions<T extends { seq: number }> {
  open: (request: ResumableStreamRequest) => Promise<{ stream: AsyncIterable<T> }>;
  from?: number | undefined;
  signal?: AbortSignal | undefined;
  maxAttempts?: number | undefined;
  retryDelay?: number | undefined;
  maxRetryDelay?: number | undefined;
  sleep?: ((ms: number, signal: AbortSignal | undefined) => Promise<void>) | undefined;
}

export class StreamEndedError extends Error {
  constructor(attempts: number, cause: unknown) {
    super(`the event stream ended ${attempts} times in a row without sending an event`, { cause });
    this.name = "StreamEndedError";
  }
}

export async function* resumableStream<T extends { seq: number }>(
  options: ResumableStreamOptions<T>,
): AsyncGenerator<T> {
  const { signal } = options;
  const sleep = options.sleep ?? abortableDelay;
  const retryDelay = options.retryDelay ?? SSE_RETRY_DELAY_MS;
  const maxRetryDelay = options.maxRetryDelay ?? SSE_MAX_RETRY_DELAY_MS;
  let from = options.from;
  let failures = 0;

  while (!isAborted(signal)) {
    const connection = new AbortController();
    const abort = () => connection.abort(signal?.reason);
    let failure: unknown;
    let received = false;

    signal?.addEventListener("abort", abort, { once: true });

    try {
      const { stream } = await options.open({
        query: from === undefined ? {} : { from },
        signal: connection.signal,
        sseMaxRetryAttempts: 1,
        onSseError: (error) => {
          failure = error;

          if (isFatalStreamError(error)) {
            throw error;
          }
        },
      });

      const events = stream[Symbol.asyncIterator]();

      try {
        for (let next = await events.next(); next.done !== true; next = await events.next()) {
          received = true;
          from = next.value.seq + 1;
          yield next.value;
        }
      } finally {
        connection.abort();
        await events.return?.();
      }
    } finally {
      signal?.removeEventListener("abort", abort);
    }

    if (isAborted(signal)) {
      return;
    }

    failures = received ? 0 : failures + 1;

    if (options.maxAttempts !== undefined && failures >= options.maxAttempts) {
      throw new StreamEndedError(failures, failure);
    }

    await sleep(Math.min(retryDelay * 2 ** Math.max(failures - 1, 0), maxRetryDelay), signal);
  }
}

function isFatalStreamError(error: unknown): error is ZodError | ApiError {
  return error instanceof ZodError || (error instanceof ApiError && error.status < 500);
}

function isAborted(signal: AbortSignal | undefined): boolean {
  return signal?.aborted === true;
}

function abortableDelay(ms: number, signal: AbortSignal | undefined): Promise<void> {
  return new Promise((resolve) => {
    const done = () => {
      clearTimeout(timer);
      signal?.removeEventListener("abort", done);
      resolve();
    };

    const timer = setTimeout(done, ms);

    signal?.addEventListener("abort", done, { once: true });
  });
}

function rejectErrorResponses(fetchFn: typeof globalThis.fetch): typeof globalThis.fetch {
  return async (input, init) => {
    const response = await fetchFn(input, init);

    if (!response.ok) {
      throw new ApiError(response.status, await errorResponse(response));
    }

    return response;
  };
}

async function errorResponse(response: Response): Promise<ErrorResponse> {
  const parsed = zErrorResponse.safeParse(
    await response
      .clone()
      .json()
      .catch(() => undefined),
  );

  if (parsed.success) {
    return parsed.data;
  }

  const text = await response.text();
  const message = text === "" ? `${response.status} ${response.statusText}`.trim() : text;

  return { code: "internal_error", message, causes: [] };
}
