import { useTheme } from "next-themes";
import { useCallback, type ComponentType, type SVGProps } from "react";
import Monitor from "~icons/lucide/monitor";
import Moon from "~icons/lucide/moon";
import Sun from "~icons/lucide/sun";

import { m, useLocale, type Locale } from "../i18n/locale";
import { Button } from "../ui/button";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const THEME_CYCLE = ["system", "light", "dark"] as const;

export type ThemeChoice = (typeof THEME_CYCLE)[number];

const THEME_META: Record<
  ThemeChoice,
  { label: (locale: Locale) => string; icon: ComponentType<SVGProps<SVGSVGElement>> }
> = {
  system: { label: (locale) => m.chrome_theme_system({}, { locale }), icon: Monitor },
  light: { label: (locale) => m.chrome_theme_light({}, { locale }), icon: Sun },
  dark: { label: (locale) => m.chrome_theme_dark({}, { locale }), icon: Moon },
};

export function nextTheme(theme: string | undefined): ThemeChoice {
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";

  return THEME_CYCLE[(THEME_CYCLE.indexOf(choice) + 1) % THEME_CYCLE.length] ?? "system";
}

export function ThemeCycleButton({ onChosen }: { onChosen?: (choice: ThemeChoice) => void }) {
  const locale = useLocale();
  const { theme, setTheme } = useTheme();
  const choice = THEME_CYCLE.find((candidate) => candidate === theme) ?? "system";
  const { label: labelOf, icon: Icon } = THEME_META[choice];
  const label = labelOf(locale);

  const cycleTheme = useCallback(() => {
    const next = nextTheme(theme);
    setTheme(next);
    onChosen?.(next);
  }, [onChosen, setTheme, theme]);

  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button variant="ghost" size="icon" aria-label={label} onClick={cycleTheme}>
            <Icon aria-hidden />
          </Button>
        }
      />
      <TooltipContent>{m.chrome_theme_tooltip({ label }, { locale })}</TooltipContent>
    </Tooltip>
  );
}
