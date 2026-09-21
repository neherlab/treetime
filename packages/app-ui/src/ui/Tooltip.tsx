import { Tooltip as BaseTooltip } from "@base-ui-components/react/tooltip";

import { cn } from "./cn";

const TooltipProvider = BaseTooltip.Provider;

const TooltipRoot = BaseTooltip.Root;

const TooltipTrigger = BaseTooltip.Trigger;

function TooltipPopup({ className, children, ...props }: BaseTooltip.Popup.Props) {
  return (
    <BaseTooltip.Portal>
      <BaseTooltip.Positioner sideOffset={6}>
        <BaseTooltip.Popup
          className={cn(
            "border-line bg-surface-3 text-2xs text-ink rounded-md border px-2 py-1 shadow-md outline-none",
            className,
          )}
          {...props}
        >
          {children}
        </BaseTooltip.Popup>
      </BaseTooltip.Positioner>
    </BaseTooltip.Portal>
  );
}

export const Tooltip = {
  Provider: TooltipProvider,
  Root: TooltipRoot,
  Trigger: TooltipTrigger,
  Popup: TooltipPopup,
};
