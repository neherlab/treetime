import { Switch as BaseSwitch } from "@base-ui/react/switch";

import { cn } from "./cn";

export interface SwitchProps {
  checked: boolean;
  onCheckedChange: (checked: boolean) => void;
  label?: string | undefined;
  id?: string | undefined;
  className?: string | undefined;
}

export function Switch({ checked, onCheckedChange, label, id, className }: SwitchProps) {
  return (
    <BaseSwitch.Root
      id={id}
      checked={checked}
      onCheckedChange={onCheckedChange}
      aria-label={label}
      className={cn(
        "bg-line-strong data-[checked]:bg-accent focus-visible:ring-accent relative inline-flex h-5 w-9 shrink-0 cursor-pointer items-center rounded-full transition-colors outline-none focus-visible:ring-2",
        className,
      )}
    >
      <BaseSwitch.Thumb className="block size-3.5 translate-x-0.5 rounded-full bg-white shadow-sm transition-transform data-[checked]:translate-x-[1.125rem]" />
    </BaseSwitch.Root>
  );
}
