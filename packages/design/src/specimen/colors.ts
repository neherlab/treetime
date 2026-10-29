import { contrastRatio, hexPair, type HexPair } from "../tokens/color";
import { CONTRAST_REQUIREMENTS, resolveColor, type BaseColorName, type ColorName } from "../tokens/tokens";

export const COLOR_GROUPS: readonly ColorGroup[] = [
  {
    title: "Surfaces",
    rows: [
      { name: "background", swatch: "bg-background", use: "Ground of the app and of the docs" },
      { name: "card", swatch: "bg-card", use: "Panels, tables, articles" },
      { name: "popover", swatch: "bg-popover", use: "Menus, popovers, dialogs" },
      { name: "muted", swatch: "bg-muted", use: "Hover, tracks, code, secondary buttons" },
      { name: "selected", swatch: "bg-selected", use: "Selected rows and items" },
    ],
  },
  {
    title: "Text and lines",
    rows: [
      { name: "foreground", swatch: "bg-foreground", use: "Text, headings, values" },
      { name: "muted-foreground", swatch: "bg-muted-foreground", use: "Descriptions, table headers, dates" },
      { name: "border", swatch: "bg-border", use: "Lines between rows, panel borders" },
      { name: "input", swatch: "bg-input", use: "Outlines of inputs and checkboxes" },
      { name: "ring", swatch: "bg-ring", use: "Keyboard focus" },
    ],
  },
  {
    title: "Actions and identity",
    rows: [
      {
        name: "primary",
        swatch: "bg-primary",
        use: "Primary buttons, links, selection; the link blue of neherlab.org",
      },
      { name: "brand", swatch: "bg-brand", use: "Page titles of sites and docs; the @ of a wordmark" },
      { name: "mark-gray", swatch: "bg-mark-gray", use: "Lines of logos" },
      { name: "mark-coral", swatch: "bg-mark-coral", use: "One stroke of a logo" },
      { name: "mark-blue", swatch: "bg-mark-blue", use: "One tip of a logo" },
    ],
  },
  {
    title: "Status",
    rows: [
      { name: "destructive", swatch: "bg-destructive", use: "Errors and destructive actions, with an icon" },
      { name: "warning", swatch: "bg-warning", use: "Warnings, with an icon" },
      { name: "success", swatch: "bg-success", use: "Finished, converged, valid" },
    ],
  },
  {
    title: "Data",
    rows: [
      { name: "data-blue", swatch: "bg-data-blue", use: "Paul Tol's bright scheme, first series" },
      { name: "data-red", swatch: "bg-data-red", use: "Second series; outliers" },
      { name: "data-green", swatch: "bg-data-green", use: "Third series" },
      { name: "data-yellow", swatch: "bg-data-yellow", use: "Fourth series; label it, below 3:1 on white" },
      { name: "data-cyan", swatch: "bg-data-cyan", use: "Fifth series; label it, below 3:1 on white" },
      { name: "data-purple", swatch: "bg-data-purple", use: "Sixth series; selection in charts" },
      { name: "data-gray", swatch: "bg-data-gray", use: "Seventh series; label it, below 3:1 on white" },
    ],
  },
];

export function colorFacts(name: ColorName): ColorFacts {
  return { hex: hexPair(resolveColor(name)), weakest: weakestContrast(name) };
}

export interface ColorGroup {
  title: string;
  rows: readonly ColorRow[];
}

export interface ColorRow {
  name: ColorName;
  swatch: string;
  use: string;
}

export interface ColorFacts {
  hex: HexPair;
  weakest: WeakestContrast | undefined;
}

export interface WeakestContrast {
  partner: BaseColorName;
  minimum: number;
  light: number;
  dark: number;
}

function weakestContrast(name: ColorName): WeakestContrast | undefined {
  const checks = CONTRAST_REQUIREMENTS.filter(
    (requirement) => requirement.foreground === name || requirement.background === name,
  ).map((requirement) => {
    const foreground = resolveColor(requirement.foreground);
    const background = resolveColor(requirement.background);

    return {
      partner: requirement.foreground === name ? requirement.background : requirement.foreground,
      minimum: requirement.minimum,
      light: contrastRatio(foreground.light, background.light),
      dark: contrastRatio(foreground.dark, background.dark),
    };
  });

  return checks.reduce<WeakestContrast | undefined>(
    (weakest, check) =>
      weakest === undefined || Math.min(check.light, check.dark) < Math.min(weakest.light, weakest.dark)
        ? check
        : weakest,
    undefined,
  );
}
