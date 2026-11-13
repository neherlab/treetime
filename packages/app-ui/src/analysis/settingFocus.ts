import type { UiDraft } from "@neherlab/app-contracts";
import { useEffect } from "react";

import { useDraftStore } from "../store/draft";
import { settingFieldId } from "./fieldIds";

const FOCUSABLE = "input, select, textarea, button, [tabindex]:not([tabindex='-1'])";

export function useSettingFocus(): void {
  const focus = useDraftStore((state) => state.focus);
  const view = useDraftStore((state) => state.draft.view);
  const search = useDraftStore((state) => state.draft.search);
  const changedOnly = useDraftStore((state) => state.draft.changed_only);

  useEffect(() => {
    if (focus === undefined) {
      return;
    }

    const target = settingElement(focus.key);
    const step = settingFocusStep({ found: target !== null, view, filtered: search !== "" || changedOnly });
    const store = useDraftStore.getState();

    if (step === "focus" && target !== null) {
      focusElement(target);
      store.clearFocus(focus.token);
    } else if (step === "open" || step === "reveal") {
      store.update({ view: "all", search: "", changed_only: false });
    }
  }, [changedOnly, focus, search, view]);
}

export function initialSettingFocus(): HTMLElement | true {
  const { focus, clearFocus } = useDraftStore.getState();
  const target = focus === undefined ? null : settingElement(focus.key);

  if (focus === undefined) {
    return true;
  }

  clearFocus(focus.token);
  target?.scrollIntoView({ block: "center" });

  return target ?? true;
}

export function settingFocusStep({ found, view, filtered }: FocusSituation): FocusStep {
  if (found) {
    return "focus";
  }

  if (view !== "all") {
    return "open";
  }

  return filtered ? "reveal" : "wait";
}

export interface FocusSituation {
  found: boolean;
  view: UiDraft["view"];
  filtered: boolean;
}

export type FocusStep = "focus" | "open" | "reveal" | "wait";

function settingElement(key: string): HTMLElement | null {
  return focusableElement(document.getElementById(settingFieldId(key, "all")));
}

function focusElement(target: HTMLElement): void {
  target.focus();
  target.scrollIntoView({ block: "center" });
}

function focusableElement(element: HTMLElement | null): HTMLElement | null {
  if (element === null) {
    return null;
  }

  return element.matches(FOCUSABLE) ? element : (element.parentElement?.querySelector<HTMLElement>(FOCUSABLE) ?? null);
}
