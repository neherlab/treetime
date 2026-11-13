import { toast } from "@neherlab/design/ui/toast";
import { useQuery } from "@tanstack/react-query";
import { useTheme } from "next-themes";
import { useEffect, useRef } from "react";

import { lostFilesMessage, reloadedFilesMessage } from "../analysis/lostRuns";
import { useApiContext } from "../api/context";
import { LoadingState } from "../components/PageShell";
import { useHost } from "../host-context";
import { reportToConsole, type PreferencesStorage } from "./storage";
import { applyPreferences, PreferencesSaver, restorePreferences } from "./sync";

const PREFERENCES_KEY = ["preferences"] as const;

export function PreferencesProvider({ storage, children }: { storage: PreferencesStorage; children: React.ReactNode }) {
  const { setTheme } = useTheme();
  const { client } = useApiContext();
  const host = useHost();

  const { data: restored, error } = useQuery({
    queryKey: PREFERENCES_KEY,
    queryFn: async () => {
      const { preferences, dropped, examples } = await restorePreferences(
        client,
        await storage.load(),
        reportToConsole,
      );

      applyPreferences(preferences, examples, setTheme);

      return { dropped };
    },
    staleTime: Infinity,
    gcTime: Infinity,
    retry: false,
  });

  const loaded = restored !== undefined;
  const dropped = restored?.dropped;
  const reported = useRef(false);

  useEffect(() => {
    if (!reported.current && dropped !== undefined && dropped.length > 0) {
      reported.current = true;
      toast.add(host === null ? reloadedFilesMessage(dropped) : lostFilesMessage(dropped));
    }
  }, [dropped, host]);

  useEffect(() => {
    if (!loaded) {
      return undefined;
    }

    const saver = new PreferencesSaver(storage, globalThis.window, reportSaveFailure);

    return () => {
      saver.dispose();
    };
  }, [loaded, storage]);

  if (error !== null) {
    throw error;
  }

  return loaded ? children : <LoadingState text="Starting TreeTime" />;
}

function reportSaveFailure(message: string): void {
  reportToConsole(`the preferences could not be saved: ${message}`);
  toast.add({ title: "The preferences could not be saved", description: message });
}
