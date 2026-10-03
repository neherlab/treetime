import { useQuery } from "@tanstack/react-query";
import { useTheme } from "next-themes";
import { useEffect } from "react";

import { LoadingState } from "../components/PageShell";
import { toast } from "../ui/toast";
import type { PreferencesStorage } from "./storage";
import { applyPreferences, PreferencesSaver } from "./sync";

const PREFERENCES_KEY = ["preferences"] as const;

export function PreferencesProvider({ storage, children }: { storage: PreferencesStorage; children: React.ReactNode }) {
  const { setTheme } = useTheme();

  const { data: loaded, error } = useQuery({
    queryKey: PREFERENCES_KEY,
    queryFn: async () => {
      applyPreferences(await storage.load(), setTheme);

      return true;
    },
    staleTime: Infinity,
    gcTime: Infinity,
    retry: false,
  });

  useEffect(() => {
    if (loaded !== true) {
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

  return loaded === true ? children : <LoadingState text="Starting TreeTime" />;
}

function reportSaveFailure(message: string): void {
  console.error(`[TreeTime] the preferences could not be saved: ${message}`);
  toast.add({ title: "The preferences could not be saved", description: message });
}
