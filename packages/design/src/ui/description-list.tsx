import { cn } from "cn";
import type { ComponentProps } from "react";

function DescriptionList({ className, ...props }: ComponentProps<"dl">) {
  return (
    <dl
      data-slot="description-list"
      className={cn("grid grid-cols-[max-content_1fr] gap-x-5 gap-y-1.5 text-sm", className)}
      {...props}
    />
  );
}

function DescriptionTerm({ className, ...props }: ComponentProps<"dt">) {
  return <dt data-slot="description-term" className={cn("text-muted-foreground", className)} {...props} />;
}

function DescriptionDetails({ className, ...props }: ComponentProps<"dd">) {
  return <dd data-slot="description-details" className={cn("text-foreground min-w-0", className)} {...props} />;
}

export { DescriptionList, DescriptionTerm, DescriptionDetails };
