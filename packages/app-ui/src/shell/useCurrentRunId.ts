import { useRouterState } from "@tanstack/react-router";

const RUN_PATH = /^\/runs\/([^/]+)/u;

export function useCurrentRunId(): string | null {
  return useRouterState({
    select: (state) => {
      const id = RUN_PATH.exec(state.location.pathname)?.[1];

      return id === undefined ? null : decodeURIComponent(id);
    },
  });
}
