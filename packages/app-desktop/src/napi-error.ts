import { errorMessage, zErrorResponse, type ErrorResponse } from "@neherlab/app-contracts";

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- the addon throws untyped errors whose message is a JSON error response, parsed here at the boundary
export function napiErrorResponse(error: unknown): ErrorResponse {
  const message = errorMessage(error);

  try {
    return zErrorResponse.parse(JSON.parse(message));
  } catch {
    return { code: "internal_error", message, causes: [] };
  }
}
