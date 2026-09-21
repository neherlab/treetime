import { Dialog as BaseDialog } from "@base-ui-components/react/dialog";

import { cn } from "./cn";

const DialogRoot = BaseDialog.Root;

const DialogTrigger = BaseDialog.Trigger;

const DialogClose = BaseDialog.Close;

const DialogTitle = BaseDialog.Title;

const DialogDescription = BaseDialog.Description;

function DialogPopup({ className, children, ...props }: BaseDialog.Popup.Props) {
  return (
    <BaseDialog.Portal>
      <BaseDialog.Backdrop className="bg-ink/25 fixed inset-0 transition-opacity data-[ending-style]:opacity-0 data-[starting-style]:opacity-0" />
      <BaseDialog.Popup
        className={cn(
          "border-line bg-surface-1 fixed top-1/2 left-1/2 w-full max-w-md -translate-x-1/2 -translate-y-1/2 rounded-lg border p-5 shadow-lg outline-none",
          className,
        )}
        {...props}
      >
        {children}
      </BaseDialog.Popup>
    </BaseDialog.Portal>
  );
}

export const Dialog = {
  Root: DialogRoot,
  Trigger: DialogTrigger,
  Close: DialogClose,
  Popup: DialogPopup,
  Title: DialogTitle,
  Description: DialogDescription,
};
