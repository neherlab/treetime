import type { RunStatus } from "@neherlab/app-contracts";
import { Ban, CircleCheck, CircleDashed, CircleX, LoaderCircle, OctagonAlert, type LucideIcon } from "lucide-react";

import { cn } from "../ui";

const STATUS_LOOK: Record<RunStatus, { label: string; icon: LucideIcon; className: string }> = {
  created: { label: "Waiting to start", icon: CircleDashed, className: "text-ink-faint" },
  running: { label: "Running", icon: LoaderCircle, className: "text-accent animate-spin" },
  ok: { label: "Finished", icon: CircleCheck, className: "text-signal-ok" },
  error: { label: "Failed", icon: CircleX, className: "text-signal-danger" },
  cancelled: { label: "Cancelled", icon: Ban, className: "text-ink-faint" },
  interrupted: { label: "Interrupted", icon: OctagonAlert, className: "text-signal-warn" },
};

export function statusLabel(status: RunStatus): string {
  return STATUS_LOOK[status].label;
}

export function StatusIcon({ status, size = 14 }: { status: RunStatus; size?: number }) {
  const look = STATUS_LOOK[status];
  const Icon = look.icon;

  return <Icon size={size} aria-label={look.label} className={cn(look.className)} />;
}
