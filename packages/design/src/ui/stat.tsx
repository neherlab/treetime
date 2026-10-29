import { cn } from "cn";
import type { ComponentProps, ReactNode } from "react";

function Stat({
  label,
  value,
  uncertainty,
  detail,
  className,
  ...props
}: ComponentProps<"div"> & { label: ReactNode; value: ReactNode; uncertainty?: ReactNode; detail?: ReactNode }) {
  return (
    <div data-slot="stat" className={cn("flex min-w-0 flex-col gap-0.5", className)} {...props}>
      <span className="text-muted-foreground text-xs">{label}</span>
      <span className="text-foreground text-xl font-medium tracking-tight tabular-nums">
        {value}
        {uncertainty !== undefined && (
          <span className="text-muted-foreground ml-1 text-sm font-normal tracking-normal">{uncertainty}</span>
        )}
      </span>
      {detail !== undefined && <span className="text-muted-foreground text-xs">{detail}</span>}
    </div>
  );
}

export { Stat };
