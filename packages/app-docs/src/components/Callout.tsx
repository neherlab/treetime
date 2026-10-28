import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import type { ComponentType, ReactNode, SVGProps } from "react";
import CircleAlert from "~icons/lucide/circle-alert";
import Info from "~icons/lucide/info";
import Lightbulb from "~icons/lucide/lightbulb";
import TriangleAlert from "~icons/lucide/triangle-alert";

const CALLOUTS = {
  info: { variant: "default", icon: Info },
  success: { variant: "success", icon: Lightbulb },
  warning: { variant: "warning", icon: TriangleAlert },
  error: { variant: "destructive", icon: CircleAlert },
} satisfies Record<string, Callout>;

export function CalloutContainer({ type, children }: { type?: string; children?: ReactNode }) {
  const callout = calloutOf(type);
  const Icon = callout.icon;

  return (
    <Alert role="note" variant={callout.variant} className="not-prose my-6">
      <Icon aria-hidden />
      {children}
    </Alert>
  );
}

export function CalloutTitle({ children }: { children?: ReactNode }) {
  return <AlertTitle>{children}</AlertTitle>;
}

export function CalloutDescription({ children }: { children?: ReactNode }) {
  return <AlertDescription className="text-foreground [&_p]:my-1">{children}</AlertDescription>;
}

interface Callout {
  variant: "default" | "success" | "warning" | "destructive";
  icon: ComponentType<SVGProps<SVGSVGElement>>;
}

function calloutOf(type: string | undefined): Callout {
  const callouts: ReadonlyMap<string, Callout> = new Map(Object.entries(CALLOUTS));

  return callouts.get(type ?? "info") ?? CALLOUTS.info;
}
