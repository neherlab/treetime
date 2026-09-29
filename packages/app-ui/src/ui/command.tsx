import { Autocomplete } from "@base-ui/react/autocomplete";
import type * as React from "react";
import SearchIcon from "~icons/lucide/search";

import { cn } from "./cn";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "./dialog";

function Command<Groups extends readonly { items: readonly unknown[] }[]>(
  props: Omit<Autocomplete.Root.Props<Groups[number]["items"][number]>, "items" | "filteredItems"> & {
    items: Groups;
    filteredItems?: Groups;
  },
) {
  return <Autocomplete.Root open inline autoHighlight="always" keepHighlight {...props} />;
}

function CommandDialog({
  title,
  description,
  children,
  className,
  finalFocus,
  ...props
}: Omit<React.ComponentProps<typeof Dialog>, "children"> & {
  title: string;
  description: string;
  className?: string;
  finalFocus?: React.ComponentProps<typeof DialogContent>["finalFocus"];
  children: React.ReactNode;
}) {
  return (
    <Dialog {...props}>
      <DialogContent
        className={cn("top-1/3 flex translate-y-0 flex-col gap-0 overflow-hidden p-0", className)}
        showCloseButton={false}
        finalFocus={finalFocus}
      >
        <DialogHeader className="sr-only">
          <DialogTitle>{title}</DialogTitle>
          <DialogDescription>{description}</DialogDescription>
        </DialogHeader>
        {children}
      </DialogContent>
    </Dialog>
  );
}

function CommandInput({ className, ...props }: Autocomplete.Input.Props) {
  return (
    <Autocomplete.InputGroup data-slot="command-input-wrapper" className="flex h-11 items-center gap-2 border-b px-3">
      <SearchIcon aria-hidden className="text-muted-foreground size-4 shrink-0" />
      <Autocomplete.Input
        data-slot="command-input"
        className={cn(
          "placeholder:text-muted-foreground h-full w-full bg-transparent text-sm outline-hidden disabled:cursor-not-allowed disabled:opacity-50",
          className,
        )}
        {...props}
      />
    </Autocomplete.InputGroup>
  );
}

function CommandList({ className, ...props }: Autocomplete.List.Props) {
  return (
    <Autocomplete.List
      data-slot="command-list"
      className={cn(
        "no-scrollbar max-h-72 scroll-py-1 overflow-x-hidden overflow-y-auto overscroll-contain p-1 outline-none",
        className,
      )}
      {...props}
    />
  );
}

function CommandEmpty({ className, ...props }: Autocomplete.Empty.Props) {
  return (
    <Autocomplete.Empty
      data-slot="command-empty"
      className={cn("text-muted-foreground py-6 text-center text-sm empty:hidden", className)}
      {...props}
    />
  );
}

function CommandGroup({ className, ...props }: Autocomplete.Group.Props) {
  return (
    <Autocomplete.Group
      data-slot="command-group"
      className={cn("text-foreground overflow-hidden not-last:mb-1", className)}
      {...props}
    />
  );
}

function CommandGroupLabel({ className, ...props }: Autocomplete.GroupLabel.Props) {
  return (
    <Autocomplete.GroupLabel
      data-slot="command-group-label"
      className={cn("text-muted-foreground px-2 py-1.5 text-xs font-bold", className)}
      {...props}
    />
  );
}

const CommandCollection = Autocomplete.Collection;

function CommandItem({ className, ...props }: Autocomplete.Item.Props) {
  return (
    <Autocomplete.Item
      data-slot="command-item"
      className={cn(
        "data-highlighted:bg-accent data-highlighted:text-accent-foreground relative flex cursor-default items-center gap-2 rounded-md px-2 py-1.5 text-sm outline-hidden select-none data-disabled:pointer-events-none data-disabled:opacity-50 [&_svg]:pointer-events-none [&_svg]:shrink-0 [&_svg:not([class*='size-'])]:size-4",
        className,
      )}
      {...props}
    />
  );
}

export {
  Command,
  CommandCollection,
  CommandDialog,
  CommandEmpty,
  CommandGroup,
  CommandGroupLabel,
  CommandInput,
  CommandItem,
  CommandList,
};
