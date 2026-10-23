import { version as getVersion } from "@neherlab/app-contracts/client";
import { useCallback } from "react";
import BookOpen from "~icons/lucide/book-open";

import { useApi } from "../api/hooks";
import { HeaderBar } from "../chrome/HeaderBar";
import { ProductIdentity } from "../chrome/ProductIdentity";
import { SearchTrigger } from "../chrome/SearchTrigger";
import { ThemeCycleButton } from "../chrome/ThemeCycleButton";
import { PALETTE_HOTKEY } from "../hotkeys";
import { m, useLocale } from "../i18n/locale";
import { usePreferencesStore } from "../store/preferences";
import { useShellStore } from "../store/shell";
import { Button } from "../ui/button";
import { Separator } from "../ui/separator";
import { SidebarTrigger } from "../ui/sidebar";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";
import { WorkspaceButton } from "./WorkspaceDialog";

const DOCS_URL = import.meta.env.VITE_TREETIME_DOCS_URL ?? "https://docs.treetime.neherlab.click";

export function SiteHeader() {
  const locale = useLocale();
  const { data: version } = useApi((context) => getVersion(context), { staleTime: Infinity });
  const setPaletteOpen = useShellStore((state) => state.setPaletteOpen);
  const storeTheme = usePreferencesStore((state) => state.setTheme);
  const openPalette = useCallback(() => setPaletteOpen(true), [setPaletteOpen]);

  return (
    <HeaderBar>
      <SidebarTrigger />
      <Separator orientation="vertical" className="data-vertical:h-4 data-vertical:self-auto" />
      <ProductIdentity logoSrc={`${import.meta.env.BASE_URL}logo-small.svg`} version={version?.version} />
      <SearchTrigger label={m.chrome_search_app({}, { locale })} hotkey={PALETTE_HOTKEY} onClick={openPalette} />
      <WorkspaceButton />
      <DocumentationButton />
      <ThemeCycleButton onChosen={storeTheme} />
    </HeaderBar>
  );
}

function DocumentationButton() {
  const locale = useLocale();
  const label = m.chrome_documentation({}, { locale });

  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button
            variant="ghost"
            size="icon"
            nativeButton={false}
            render={<a href={DOCS_URL} target="_blank" rel="noreferrer" aria-label={label} />}
          >
            <BookOpen aria-hidden />
          </Button>
        }
      />
      <TooltipContent>{label}</TooltipContent>
    </Tooltip>
  );
}
