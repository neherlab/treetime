import { cn } from "cn";
import {
  Group,
  type GroupProps,
  Panel,
  type PanelProps,
  type PanelSize,
  Separator,
  type SeparatorProps,
} from "react-resizable-panels";

function ResizablePanelGroup({ className, ...props }: GroupProps) {
  return (
    <Group
      data-slot="resizable-panel-group"
      className={cn("flex h-full w-full aria-[orientation=vertical]:flex-col", className)}
      {...props}
    />
  );
}

function ResizablePanel({ ...props }: PanelProps) {
  return <Panel data-slot="resizable-panel" {...props} />;
}

function ResizableHandle({ className, ...props }: SeparatorProps) {
  return (
    <Separator
      data-slot="resizable-handle"
      className={cn(
        "bg-border ring-offset-background focus-visible:ring-ring data-[separator=active]:bg-primary data-[separator=hover]:bg-primary relative flex w-px items-center justify-center transition-colors after:absolute after:inset-y-0 after:left-1/2 after:w-1 after:-translate-x-1/2 focus-visible:ring-1 focus-visible:outline-hidden aria-[orientation=horizontal]:h-px aria-[orientation=horizontal]:w-full aria-[orientation=horizontal]:after:left-0 aria-[orientation=horizontal]:after:h-1 aria-[orientation=horizontal]:after:w-full aria-[orientation=horizontal]:after:translate-x-0 aria-[orientation=horizontal]:after:-translate-y-1/2",
        className,
      )}
      {...props}
    />
  );
}

export { type PanelSize, ResizableHandle, ResizablePanel, ResizablePanelGroup };
