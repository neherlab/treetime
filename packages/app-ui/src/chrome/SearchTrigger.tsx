import { Button } from "@neherlab/design/ui/button";
import { Kbd } from "@neherlab/design/ui/kbd";
import { formatForDisplay, type Hotkey } from "@tanstack/react-hotkeys";
import Search from "~icons/lucide/search";

export function SearchTrigger({ label, hotkey, onClick }: { label: string; hotkey: Hotkey; onClick: () => void }) {
  return (
    <Button
      variant="outline"
      aria-label={label}
      onClick={onClick}
      className="text-muted-foreground ml-auto justify-start sm:w-72"
    >
      <Search aria-hidden />
      <span className="hidden flex-1 text-left sm:inline">{label}</span>
      <Kbd className="hidden sm:inline-flex">{formatForDisplay(hotkey)}</Kbd>
    </Button>
  );
}
