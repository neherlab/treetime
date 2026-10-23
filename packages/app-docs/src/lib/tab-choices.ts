import { create } from "zustand";
import { persist } from "zustand/middleware";

export const useTabChoices = create<TabChoices>()(
  persist(
    (set) => ({
      choices: {},
      choose: (syncKey, label) => {
        set((state) => ({ choices: { ...state.choices, [syncKey]: label } }));
      },
    }),
    { name: "treetime-docs-tabs", skipHydration: true },
  ),
);

interface TabChoices {
  choices: Record<string, string>;
  choose: (syncKey: string, label: string) => void;
}
