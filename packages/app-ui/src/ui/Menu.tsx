import { Menu as BaseMenu } from "@base-ui-components/react/menu";

import { cn } from "./cn";

const MenuRoot = BaseMenu.Root;

const MenuTrigger = BaseMenu.Trigger;

const MenuGroup = BaseMenu.Group;

function MenuPopup({ className, children, ...props }: BaseMenu.Popup.Props) {
  return (
    <BaseMenu.Portal>
      <BaseMenu.Positioner sideOffset={6} className="outline-none">
        <BaseMenu.Popup
          className={cn(
            "border-line bg-surface-1 text-ink min-w-40 rounded-md border p-1 text-sm shadow-lg outline-none",
            className,
          )}
          {...props}
        >
          {children}
        </BaseMenu.Popup>
      </BaseMenu.Positioner>
    </BaseMenu.Portal>
  );
}

function MenuItem({ className, ...props }: BaseMenu.Item.Props) {
  return (
    <BaseMenu.Item
      className={cn(
        "data-[highlighted]:bg-accent-subtle data-[highlighted]:text-accent flex cursor-default items-center gap-2 rounded-sm px-2 py-1.5 outline-none select-none",
        className,
      )}
      {...props}
    />
  );
}

export const Menu = {
  Root: MenuRoot,
  Trigger: MenuTrigger,
  Group: MenuGroup,
  Popup: MenuPopup,
  Item: MenuItem,
};
