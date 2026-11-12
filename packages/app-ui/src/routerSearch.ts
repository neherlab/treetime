import { retainSearchParams } from "@tanstack/react-router";

export const rootSearch = {
  search: { middlewares: [retainSearchParams(true)] },
};
