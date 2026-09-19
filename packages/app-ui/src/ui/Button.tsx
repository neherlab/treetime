import { Button as BaseButton } from "@base-ui-components/react/button";
import { cva, type VariantProps } from "class-variance-authority";

import { cn } from "./cn";

const buttonVariants = cva(
  "focus-visible:ring-accent focus-visible:ring-offset-surface-0 inline-flex shrink-0 items-center justify-center gap-1.5 rounded-md font-medium whitespace-nowrap transition-colors outline-none focus-visible:ring-2 focus-visible:ring-offset-1 disabled:pointer-events-none disabled:opacity-50",
  {
    variants: {
      variant: {
        solid: "bg-accent text-accent-fg hover:bg-accent-hover",
        outline: "border-line-strong text-ink hover:bg-surface-2 border",
        ghost: "text-ink-muted hover:bg-surface-2 hover:text-ink",
        subtle: "bg-surface-2 text-ink hover:bg-surface-3",
        danger: "bg-signal-danger text-accent-fg hover:opacity-90",
      },
      size: {
        sm: "text-2xs h-7 px-2.5",
        md: "h-9 px-4 text-sm",
        icon: "size-8 p-0",
      },
    },
    defaultVariants: {
      variant: "solid",
      size: "md",
    },
  },
);

export type ButtonProps = BaseButton.Props & VariantProps<typeof buttonVariants> & { className?: string };

export function Button({ className, variant, size, ...props }: ButtonProps) {
  return <BaseButton className={cn(buttonVariants({ variant, size }), className)} {...props} />;
}

export { buttonVariants };
