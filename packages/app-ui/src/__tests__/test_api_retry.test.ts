import { ApiError } from "@neherlab/app-contracts/client";

import { retryQuery } from "../api/retry";

const NOT_FOUND = new ApiError(404, { code: "not_found", message: "no run with id r1", causes: [] });

const INTERNAL = new ApiError(500, { code: "internal_error", message: "the back end failed", causes: [] });

describe("query retry", () => {
  test.each([
    { name: "not found on the first failure", error: NOT_FOUND, failures: 0, expected: false },
    { name: "internal error on the first failure", error: INTERNAL, failures: 0, expected: true },
    {
      name: "network error on the first failure",
      error: new TypeError("failed to fetch"),
      failures: 0,
      expected: true,
    },
    { name: "internal error after one retry", error: INTERNAL, failures: 1, expected: false },
  ])("$name: $expected", ({ error, failures, expected }) => {
    expect(retryQuery(failures, error)).toBe(expected);
  });
});
