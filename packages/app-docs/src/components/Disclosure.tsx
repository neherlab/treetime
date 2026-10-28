import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import type { ReactNode } from "react";
import ChevronRight from "~icons/lucide/chevron-right";

export function Disclosure({ title, children }: { title: string; children?: ReactNode }) {
  return (
    <Collapsible className="not-prose bg-card my-6 rounded-lg border">
      <CollapsibleTrigger className="group/trigger flex w-full items-center gap-2 px-4 py-3 text-left text-sm font-bold">
        <ChevronRight aria-hidden className="size-4 transition-transform group-data-panel-open/trigger:rotate-90" />
        {title}
      </CollapsibleTrigger>
      <CollapsibleContent hiddenUntilFound className="prose px-4 pb-4">
        {children}
      </CollapsibleContent>
    </Collapsible>
  );
}
