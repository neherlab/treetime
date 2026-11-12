import { useEffect } from "react";

import { useDraftStore } from "../store/draft";
import { settingFieldId } from "./fieldIds";

const FOCUSABLE = "input, select, textarea, button, [tabindex]:not([tabindex='-1'])";

export function useSettingFocus(): void {
  const focus = useDraftStore((state) => state.focus);
  const view = useDraftStore((state) => state.draft.view);

  useEffect(() => {
    if (focus === undefined) {
      return;
    }

    const target = focusableElement(document.getElementById(settingFieldId(focus.key)));

    if (target !== null) {
      target.focus();
      target.scrollIntoView({ block: "center" });
      useDraftStore.getState().clearFocus(focus.token);
    } else if (view !== "all") {
      useDraftStore.getState().update({ view: "all", search: "", changed_only: false });
    }
  }, [focus, view]);
}

function focusableElement(element: HTMLElement | null): HTMLElement | null {
  if (element === null) {
    return null;
  }

  return element.matches(FOCUSABLE) ? element : (element.parentElement?.querySelector<HTMLElement>(FOCUSABLE) ?? null);
}
