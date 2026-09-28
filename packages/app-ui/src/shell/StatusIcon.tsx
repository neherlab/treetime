import type { RunStatus } from "@neherlab/app-contracts";
import { Ban, CircleCheck, CircleDashed, CircleX, LoaderCircle, OctagonAlert, type LucideIcon } from "lucide-react";

import { cn } from "../ui/cn";

const STATUS_LOOK: Record<RunStatus, { label: string; icon: LucideIcon; className: string }> = {
  created: { label: "Waiting to start", icon: CircleDashed, className: "text-muted-foreground" },
  running: { label: "Running", icon: LoaderCircle, className: "text-primary animate-spin" },
  ok: { label: "Finished", icon: CircleCheck, className: "text-success" },
  error: { label: "Failed", icon: CircleX, className: "text-destructive" },
  cancelled: { label: "Cancelled", icon: Ban, className: "text-muted-foreground" },
  interrupted: { label: "Interrupted", icon: OctagonAlert, className: "text-warning" },
};

export function statusLabel(status: RunStatus): string {
  return STATUS_LOOK[status].label;
}

export function StatusIcon({ status, className }: { status: RunStatus; className?: string }) {
  const look = STATUS_LOOK[status];
  const Icon = look.icon;

  return <Icon aria-label={look.label} className={cn("size-3.5 shrink-0", look.className, className)} />;
}
