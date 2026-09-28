import type * as React from "react";
import * as RechartsPrimitive from "recharts";
import type { TooltipPayloadEntry, TooltipValueType } from "recharts";

import { cn } from "./cn";

const INITIAL_DIMENSION = { width: 320, height: 200 } as const;

type TooltipNameType = number | string;

function ChartContainer({
  className,
  children,
  initialDimension = INITIAL_DIMENSION,
  ...props
}: React.ComponentProps<"div"> & {
  children: React.ComponentProps<typeof RechartsPrimitive.ResponsiveContainer>["children"];
  initialDimension?: {
    width: number;
    height: number;
  };
}) {
  return (
    <div
      data-slot="chart"
      className={cn(
        "[&_.recharts-cartesian-axis-tick_text]:fill-muted-foreground [&_.recharts-curve.recharts-tooltip-cursor]:stroke-border [&_.recharts-rectangle.recharts-tooltip-cursor]:fill-muted flex aspect-video justify-center text-xs [&_.recharts-dot[stroke='#fff']]:stroke-transparent [&_.recharts-layer]:outline-hidden [&_.recharts-surface]:outline-hidden",
        className,
      )}
      {...props}
    >
      <RechartsPrimitive.ResponsiveContainer initialDimension={initialDimension}>
        {children}
      </RechartsPrimitive.ResponsiveContainer>
    </div>
  );
}

const ChartTooltip = RechartsPrimitive.Tooltip;

function ChartTooltipFrame({ className, ...props }: React.ComponentProps<"div">) {
  return (
    <div
      data-slot="chart-tooltip"
      className={cn(
        "bg-popover text-popover-foreground grid min-w-32 items-start gap-1 rounded-lg border px-2.5 py-1.5 text-xs shadow-xl",
        className,
      )}
      {...props}
    />
  );
}

function ChartTooltipContent({
  active,
  payload,
}: {
  active?: boolean | undefined;
  payload?: readonly TooltipPayloadEntry<TooltipValueType, TooltipNameType>[] | undefined;
}) {
  if (active !== true || payload === undefined || payload.length === 0) {
    return null;
  }

  return (
    <ChartTooltipFrame>
      {payload
        .filter((item) => item.type !== "none")
        .map((item) => (
          <div key={item.graphicalItemId} className="flex w-full items-center gap-2">
            <svg viewBox="0 0 10 10" aria-hidden className="size-2.5 shrink-0">
              <rect width="10" height="10" rx="2" fill={item.color} />
            </svg>
            <div className="flex flex-1 items-center justify-between gap-2 leading-none">
              <span className="text-muted-foreground">{item.name}</span>
              {item.value !== undefined && (
                <span className="text-foreground font-mono font-medium tabular-nums">{tooltipValue(item.value)}</span>
              )}
            </div>
          </div>
        ))}
    </ChartTooltipFrame>
  );
}

function tooltipValue(value: TooltipValueType): string {
  return value.toLocaleString();
}

export { ChartContainer, ChartTooltip, ChartTooltipContent, ChartTooltipFrame };
