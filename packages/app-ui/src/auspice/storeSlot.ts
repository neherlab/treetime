import type { AuspiceDocument } from "@neherlab/app-contracts";
import { createContext, use } from "react";

import type { AuspiceStore } from "./store";

export class AuspiceStoreSlot {
  #entry: { document: AuspiceDocument; store: AuspiceStore } | undefined;

  store(document: AuspiceDocument, create: () => AuspiceStore): AuspiceStore {
    if (this.#entry?.document !== document) {
      this.#entry = { document, store: create() };
    }

    return this.#entry.store;
  }
}

export const AuspiceStoreSlotContext = createContext<AuspiceStoreSlot | undefined>(undefined);

export function useAuspiceStoreSlot(): AuspiceStoreSlot | undefined {
  return use(AuspiceStoreSlotContext);
}
