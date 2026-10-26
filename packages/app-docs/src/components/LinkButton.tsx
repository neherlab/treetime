import { buttonVariants } from "@neherlab/app-ui/ui/button";
import type { ReactNode } from "react";

export function LinkButton({ href, variant, size, children }: LinkButtonProps) {
  return (
    <a href={href} className={buttonVariants({ variant, size })}>
      {children}
    </a>
  );
}

type LinkButtonProps = { href: string; children?: ReactNode } & NonNullable<Parameters<typeof buttonVariants>[0]>;
