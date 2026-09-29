import type { RunStatus } from "@neherlab/app-contracts";
import type { ComponentType, SVGProps } from "react";
import Ban from "~icons/lucide/ban";
import CircleCheck from "~icons/lucide/circle-check";
import CircleDashed from "~icons/lucide/circle-dashed";
import CircleX from "~icons/lucide/circle-x";
import LoaderCircle from "~icons/lucide/loader-circle";
import OctagonAlert from "~icons/lucide/octagon-alert";

import { cn } from "../ui/cn";

const STATUS_LOOK: Record<
  RunStatus,
  { label: string; icon: ComponentType<SVGProps<SVGSVGElement>>; className: string }
> = {
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
