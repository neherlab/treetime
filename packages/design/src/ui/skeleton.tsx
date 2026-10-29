import { cn } from "cn";
import type { ComponentProps } from "react";

function Skeleton({ className, ...props }: ComponentProps<"div">) {
  return <div data-slot="skeleton" className={cn("bg-muted animate-pulse rounded-md", className)} {...props} />;
}

export { Skeleton };
