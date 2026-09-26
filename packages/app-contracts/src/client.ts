import { ZodError } from "zod";

import { type Client, type Config, createClient, createConfig, type RequestOptions } from "./generated/client";
import type { ErrorResponse } from "./generated/types.gen";
import { zErrorResponse } from "./generated/zod.gen";

export * from "./generated/sdk.gen";

export type { Client as ApiClient } from "./generated/client";

export const SSE_MAX_RETRY_ATTEMPTS = 8;

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
      if (error instanceof ZodError || (error instanceof ApiError && error.status < 500)) {
        throw error;
      }
    },
    sseMaxRetryAttempts: SSE_MAX_RETRY_ATTEMPTS,
  };

  return createClient(config);
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
