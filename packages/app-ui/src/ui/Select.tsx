import { Select as BaseSelect } from "@base-ui-components/react/select";
import { ChevronsUpDown } from "lucide-react";

import { cn } from "./cn";

const SelectRoot = BaseSelect.Root;
const SelectValue = BaseSelect.Value;
const SelectGroup = BaseSelect.Group;

function SelectTrigger({ className, children, ...props }: BaseSelect.Trigger.Props) {
  return (
    <BaseSelect.Trigger
      className={cn(
        "border-line bg-surface-0 text-ink hover:border-line-strong focus-visible:ring-accent data-[popup-open]:border-accent flex h-9 w-full items-center justify-between gap-2 rounded-md border px-2.5 text-sm outline-none focus-visible:ring-2",
        className,
      )}
      {...props}
    >
      {children}
      <BaseSelect.Icon className="text-ink-faint">
        <ChevronsUpDown size={14} />
      </BaseSelect.Icon>
    </BaseSelect.Trigger>
  );
}

function SelectPopup({ className, children, ...props }: BaseSelect.Popup.Props) {
  return (
    <BaseSelect.Portal>
      <BaseSelect.Positioner sideOffset={6} className="outline-none">
        <BaseSelect.Popup
          className={cn(
            "border-line bg-surface-1 text-ink max-h-72 min-w-[var(--anchor-width)] overflow-y-auto rounded-md border p-1 text-sm shadow-lg outline-none",
            className,
          )}
          {...props}
        >
          {children}
        </BaseSelect.Popup>
      </BaseSelect.Positioner>
    </BaseSelect.Portal>
  );
}

function SelectItem({ className, ...props }: BaseSelect.Item.Props) {
  return (
    <BaseSelect.Item
      className={cn(
        "data-[highlighted]:bg-accent-subtle data-[highlighted]:text-accent flex cursor-default items-center rounded-sm px-2 py-1.5 outline-none select-none data-[selected]:font-medium",
        className,
      )}
      {...props}
    />
  );
}

export const Select = {
  Root: SelectRoot,
  Trigger: SelectTrigger,
  Value: SelectValue,
  Group: SelectGroup,
  Popup: SelectPopup,
  Item: SelectItem,
};
