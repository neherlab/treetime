import { cn } from "cn";
import type { ComponentProps } from "react";

function Wordmark({ name, place, className, ...props }: ComponentProps<"span"> & { name: string; place?: string }) {
  return (
    <span
      data-slot="wordmark"
      className={cn("text-foreground font-mono font-medium tracking-tight whitespace-nowrap", className)}
      {...props}
    >
      {name}
      {place !== undefined && (
        <>
          <span className="text-brand font-normal">@</span>
          <span className="text-muted-foreground font-normal">{place}</span>
        </>
      )}
    </span>
  );
}

export { Wordmark };
