import { REMOVE_NOTIFICATION } from "auspice/src/actions/types";
import { useEffect } from "react";

import { useToastManager } from "../ui/toast";
import type { AuspiceNotification } from "./state";
import type { AuspiceStore } from "./store";

export function useNotificationToasts(store: AuspiceStore): void {
  const toasts = useToastManager();

  useEffect(() => {
    let moving = false;

    const move = () => {
      if (moving) {
        return;
      }

      moving = true;

      for (const notification of store.getState().notifications.stack) {
        toasts.add(notificationToast(notification));
        store.dispatch({ type: REMOVE_NOTIFICATION, id: notification.id });
      }

      moving = false;
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
