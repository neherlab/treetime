import { REMOVE_NOTIFICATION } from "auspice/src/actions/types";
import { useEffect } from "react";

import { useToastManager } from "../ui/toast";
import type { AuspiceNotification } from "./state";
import type { AuspiceStore } from "./store";

export function useNotificationToasts(store: AuspiceStore): void {
  const toasts = useToastManager();

  useEffect(() => {
    const move = () => {
      const [first] = store.getState().notifications.stack;

      if (first === undefined) {
        return;
      }

      toasts.add(notificationToast(first));
      store.dispatch({ type: REMOVE_NOTIFICATION, id: first.id });
    };

    move();

    return store.subscribe(move);
  }, [store, toasts]);
}

function notificationToast({ message, details, notificationType }: AuspiceNotification): NotificationToast {
  const toast: NotificationToast = { title: message, type: notificationType };

  if (details !== "") {
    toast.description = details;
  }

  if (notificationType === "error") {
    toast.timeout = 0;
  }

  return toast;
}

interface NotificationToast {
  title: string;
  type: AuspiceNotification["notificationType"];
  description?: string;
  timeout?: number;
}
