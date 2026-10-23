import { Tabs, TabsContent, TabsList, TabsTrigger } from "@neherlab/app-ui/ui/tabs";
import { useCallback, useEffect, type ReactNode } from "react";

import { useTabChoices } from "../lib/tab-choices";

export function SyncedTabs({ syncKey, labels, ...slots }: SyncedTabsProps) {
  const chosen = useTabChoices((state) => state.choices[syncKey]);
  const choose = useTabChoices((state) => state.choose);
  useEffect(() => {
    void useTabChoices.persist.rehydrate();
  }, []);

  const value = chosen !== undefined && labels.includes(chosen) ? chosen : labels[0];

  const onValueChange = useCallback(
    (next: string) => {
      const label = labels.find((candidate) => candidate === next);

      if (label !== undefined) {
        choose(syncKey, label);
      }
    },
    [choose, labels, syncKey],
  );

  return (
    <Tabs value={value} onValueChange={onValueChange} className="not-prose my-6">
      <TabsList>
        {labels.map((label) => (
          <TabsTrigger key={label} value={label}>
            {label}
          </TabsTrigger>
        ))}
      </TabsList>
      {labels.map((label, index) => (
        <TabsContent key={label} value={label} className="prose pt-2">
          {slots[`tab-${index}`]}
        </TabsContent>
      ))}
    </Tabs>
  );
}

type SyncedTabsProps = { syncKey: string; labels: string[] } & Partial<Record<`tab-${number}`, ReactNode>>;
