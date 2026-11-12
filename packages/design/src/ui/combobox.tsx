import { Autocomplete } from "@base-ui/react/autocomplete";
import { cn } from "cn";
import SearchIcon from "~icons/lucide/search";

const Combobox = Autocomplete.Root;

function ComboboxInput({ className, ...props }: Autocomplete.Input.Props) {
  return (
    <Autocomplete.InputGroup
      data-slot="combobox-input-group"
      className={cn(
        "border-input focus-within:border-ring focus-within:ring-ring/50 flex h-8 items-center gap-2 rounded-md border bg-transparent px-2.5 shadow-xs focus-within:ring-[3px]",
        className,
      )}
    >
      <SearchIcon aria-hidden className="text-muted-foreground size-4 shrink-0" />
      <Autocomplete.Input
        data-slot="combobox-input"
        className="placeholder:text-muted-foreground h-full w-full min-w-0 bg-transparent text-sm outline-hidden"
        {...props}
      />
    </Autocomplete.InputGroup>
  );
}

function ComboboxPopup({ className, children, ...props }: Autocomplete.Popup.Props) {
  return (
    <Autocomplete.Portal>
      <Autocomplete.Positioner sideOffset={4} align="start" className="isolate z-50">
        <Autocomplete.Popup
          data-slot="combobox-popup"
          className={cn(
            "bg-popover text-popover-foreground ring-foreground/10 data-open:animate-in data-open:fade-in-0 data-closed:animate-out data-closed:fade-out-0 w-(--anchor-width) min-w-56 overflow-hidden rounded-md shadow-md ring-1 duration-100",
            className,
          )}
          {...props}
        >
          {children}
        </Autocomplete.Popup>
      </Autocomplete.Positioner>
    </Autocomplete.Portal>
  );
}

function ComboboxList({ className, ...props }: Autocomplete.List.Props) {
  return (
    <Autocomplete.List
      data-slot="combobox-list"
      className={cn("max-h-72 scroll-py-1 overflow-y-auto overscroll-contain p-1 outline-none empty:p-0", className)}
      {...props}
    />
  );
}

function ComboboxItem({ className, ...props }: Autocomplete.Item.Props) {
  return (
    <Autocomplete.Item
      data-slot="combobox-item"
      className={cn(
        "data-highlighted:bg-accent data-highlighted:text-accent-foreground flex cursor-default items-center rounded-sm px-2 py-1.5 text-sm outline-hidden select-none",
        className,
      )}
      {...props}
    />
  );
}

function ComboboxEmpty({ className, ...props }: Autocomplete.Empty.Props) {
  return (
    <Autocomplete.Empty
      data-slot="combobox-empty"
      className={cn("text-muted-foreground px-2 py-3 text-sm empty:hidden", className)}
      {...props}
    />
  );
}

export { Combobox, ComboboxEmpty, ComboboxInput, ComboboxItem, ComboboxList, ComboboxPopup };
