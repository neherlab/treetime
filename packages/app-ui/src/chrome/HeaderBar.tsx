import { cn } from "cn";
import type { ComponentProps } from "react";

export function HeaderBar({ className, ...props }: ComponentProps<"header">) {
  return (
    <header
      className={cn("bg-background z-20 flex h-(--header-height) shrink-0 items-center gap-2 border-b px-3", className)}
      {...props}
    />
  );
}
