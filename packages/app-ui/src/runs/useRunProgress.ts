import type { TreeTimeBridge } from "@neherlab/app-contracts";
import { useQueryClient, type QueryClient } from "@tanstack/react-query";
import { useEffect, useState } from "react";

import { useBridge } from "../BridgeContext";
import { RUNS_KEY } from "../queries";
import { EMPTY_PROGRESS, foldRunEvents, type RunEvent, type RunProgress } from "../results/progress";

interface RunProgressState {
  progress: RunProgress;
  failure: string | undefined;
}

interface KeyedState {
  id: string;
  value: RunProgressState;
}

const INITIAL_STATE: RunProgressState = { progress: EMPTY_PROGRESS, failure: undefined };

export function useRunProgress(id: string): RunProgressState {
  const bridge = useBridge();
  const queryClient = useQueryClient();
  const [state, setState] = useState<KeyedState>({ id, value: INITIAL_STATE });

  useEffect(() => {
    const controller = new AbortController();

    void followEvents({ bridge, queryClient, id, signal: controller.signal, update: setState });

    return () => controller.abort();
  }, [bridge, id, queryClient]);

  return state.id === id ? state.value : INITIAL_STATE;
}

async function followEvents({
  bridge,
  queryClient,
  id,
  signal,
  update,
}: {
  bridge: TreeTimeBridge;
  queryClient: QueryClient;
  id: string;
  signal: AbortSignal;
  update: (change: (current: KeyedState) => KeyedState) => void;
}): Promise<void> {
  let pending: RunEvent[] = [];
  let frame: number | undefined;

  function flush() {
    frame = undefined;
    const batch = pending;
    pending = [];
    update((current) => {
      const value = current.id === id ? current.value : INITIAL_STATE;

      return { id, value: { ...value, progress: foldRunEvents(value.progress, batch) } };
    });
  }

  signal.addEventListener("abort", () => {
    if (frame !== undefined) {
      cancelAnimationFrame(frame);
    }
  });

  try {
    await bridge.followRun(id, {
      signal,
      onEvent: (event) => {
        pending.push(event);
        frame ??= requestAnimationFrame(flush);

        if (event.type === "terminal") {
          void queryClient.invalidateQueries({ queryKey: RUNS_KEY });
        }
      },
    });
  } catch (error) {
    if (!signal.aborted) {
      const failure = error instanceof Error ? error.message : String(error);

      update((current) => ({ id, value: { ...(current.id === id ? current.value : INITIAL_STATE), failure } }));
    }
  }
}
