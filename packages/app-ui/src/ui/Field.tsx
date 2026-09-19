import { Field as BaseField } from "@base-ui-components/react/field";

import { cn } from "./cn";

const FieldRoot = BaseField.Root;

function FieldLabel({ className, ...props }: BaseField.Label.Props) {
  return <BaseField.Label className={cn("text-2xs text-ink-muted font-medium", className)} {...props} />;
}

function FieldControl({ className, ...props }: BaseField.Control.Props) {
  return (
    <BaseField.Control
      className={cn(
        "border-line bg-surface-0 text-ink placeholder:text-ink-faint hover:border-line-strong focus-visible:ring-accent data-[invalid]:border-signal-danger h-9 w-full rounded-md border px-2.5 text-sm outline-none focus-visible:ring-2",
        className,
      )}
      {...props}
    />
  );
}

function FieldDescription({ className, ...props }: BaseField.Description.Props) {
  return <BaseField.Description className={cn("text-2xs text-ink-faint", className)} {...props} />;
}

function FieldError({ className, ...props }: BaseField.Error.Props) {
  return <BaseField.Error className={cn("text-2xs text-signal-danger", className)} {...props} />;
}

export const Field = {
  Root: FieldRoot,
  Label: FieldLabel,
  Control: FieldControl,
  Description: FieldDescription,
  Error: FieldError,
};
