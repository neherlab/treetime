import { useCallback } from "react";

import { Toast } from "../ui";

export function useCopy(): (text: string, message: string) => void {
  const toasts = Toast.useToastManager();

  return useCallback(
    (text: string, message: string) => {
      navigator.clipboard.writeText(text).then(
        () => toasts.add({ title: message }),
        () => toasts.add({ title: "The browser blocked clipboard access" }),
      );
    },
    [toasts],
  );
}
