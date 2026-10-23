import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { Button } from "@neherlab/app-ui/ui/button";
import { Sheet, SheetContent, SheetHeader, SheetTitle } from "@neherlab/app-ui/ui/sheet";
import { Tooltip, TooltipContent, TooltipTrigger } from "@neherlab/app-ui/ui/tooltip";
import { useCallback, useState } from "react";
import Menu from "~icons/lucide/menu";

import type { NavNode } from "../lib/nav";
import { NavTree } from "./NavTree";

export function MobileNav({ nodes, currentUrl, locale }: { nodes: NavNode[]; currentUrl: string; locale: Locale }) {
  const [open, setOpen] = useState(false);
  const show = useCallback(() => setOpen(true), []);
  const label = m.docs_open_navigation({}, { locale });

  return (
    <LocaleProvider locale={locale}>
      <Tooltip>
        <TooltipTrigger
          render={
            <Button variant="ghost" size="icon" aria-label={label} onClick={show} className="md:hidden">
              <Menu aria-hidden />
            </Button>
          }
        />
        <TooltipContent>{label}</TooltipContent>
      </Tooltip>
      <Sheet open={open} onOpenChange={setOpen}>
        <SheetContent side="left" className="overflow-y-auto">
          <SheetHeader>
            <SheetTitle>{m.docs_navigation({}, { locale })}</SheetTitle>
          </SheetHeader>
          <div className="px-2 pb-4">
            <NavTree nodes={nodes} currentUrl={currentUrl} locale={locale} />
          </div>
        </SheetContent>
      </Sheet>
    </LocaleProvider>
  );
}
