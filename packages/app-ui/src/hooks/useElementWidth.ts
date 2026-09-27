import { useCallback, useSyncExternalStore } from "react";

export function useElementWidth(element: HTMLElement | null): number {
  const subscribe = useCallback(
    (notify: () => void) => {
      if (element === null) {
        return () => undefined;
      }

      const observer = new ResizeObserver(notify);
      observer.observe(element);

      return () => observer.disconnect();
    },
    [element],
  );

  const read = useCallback(() => Math.floor(element?.clientWidth ?? 0), [element]);

  return useSyncExternalStore(subscribe, read, read);
}
