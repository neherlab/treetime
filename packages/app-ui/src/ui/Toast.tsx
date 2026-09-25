import { Toast as BaseToast } from "@base-ui/react/toast";

function ToastList() {
  const { toasts } = BaseToast.useToastManager();

  return toasts.map((toast) => (
    <BaseToast.Root
      key={toast.id}
      toast={toast}
      className="bg-ink text-surface-1 flex max-w-md items-center gap-3 rounded-lg px-3.5 py-2.5 text-sm shadow-lg"
    >
      <BaseToast.Content className="flex-1">
        <BaseToast.Title className="font-bold" />
        <BaseToast.Description />
      </BaseToast.Content>
      <BaseToast.Action className="rounded-sm border border-current/40 px-2 py-0.5 text-xs font-bold" />
      <BaseToast.Close aria-label="Dismiss" className="text-xs opacity-70 hover:opacity-100">
        Close
      </BaseToast.Close>
    </BaseToast.Root>
  ));
}

function ToastViewport() {
  return (
    <BaseToast.Portal>
      <BaseToast.Viewport className="fixed right-4 bottom-4 z-50 flex flex-col gap-2 outline-none">
        <ToastList />
      </BaseToast.Viewport>
    </BaseToast.Portal>
  );
}

export const Toast = {
  Provider: BaseToast.Provider,
  Viewport: ToastViewport,
  useToastManager: BaseToast.useToastManager,
};
