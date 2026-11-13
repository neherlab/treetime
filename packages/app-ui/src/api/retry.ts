import { ApiError } from "@neherlab/app-contracts/client";

const QUERY_RETRIES = 1;

export function retryQuery(failureCount: number, error: Error): boolean {
  return failureCount < QUERY_RETRIES && !isNotFound(error);
}

export function isNotFound(error: Error): boolean {
  return error instanceof ApiError && error.response.code === "not_found";
}
