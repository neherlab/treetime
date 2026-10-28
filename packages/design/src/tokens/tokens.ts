const BLUE_HUE = 259.5;

export const BASE_COLORS = {
  background: pair(gray(0.985), gray(0.2)),
  foreground: pair(gray(0.21), gray(0.93)),
  card: pair(gray(1), gray(0.23)),
  popover: pair(gray(1), gray(0.26)),
  primary: pair(blue(0.53, 0.17), blue(0.55, 0.16)),
  "primary-foreground": pair(gray(1), gray(1)),
  muted: pair(gray(0.962), gray(0.28)),
  "muted-foreground": pair(gray(0.505), gray(0.72)),
  border: pair(gray(0.912), gray(0.315)),
  input: pair(gray(0.645), gray(0.54)),
  ring: pair(blue(0.56, 0.17), blue(0.72, 0.13)),
  link: pair(blue(0.51, 0.17), blue(0.76, 0.12)),
  selected: pair(blue(0.955, 0.021), blue(0.3, 0.07)),
  brand: pair(oklch(0.42, 0.155, 27), oklch(0.68, 0.16, 25)),
  destructive: pair(oklch(0.53, 0.19, 33), oklch(0.72, 0.15, 30)),
  "destructive-subtle": pair(oklch(0.965, 0.016, 30), oklch(0.28, 0.06, 30)),
  warning: pair(oklch(0.53, 0.115, 65), oklch(0.82, 0.13, 80)),
  "warning-subtle": pair(oklch(0.975, 0.025, 85), oklch(0.29, 0.05, 75)),
  success: pair(oklch(0.5, 0.12, 150), oklch(0.77, 0.14, 152)),
  "success-subtle": pair(oklch(0.97, 0.03, 150), oklch(0.28, 0.05, 152)),
  "mark-gray": pair(gray(0.573), gray(0.68)),
  "mark-coral": pair(oklch(0.622, 0.141, 22), oklch(0.7, 0.13, 22)),
  "mark-blue": pair(oklch(0.589, 0.087, 228), oklch(0.7, 0.08, 228)),
  "data-blue": pair(oklch(0.556, 0.097, 250), oklch(0.68, 0.1, 250)),
  "data-red": pair(oklch(0.683, 0.168, 15.1), oklch(0.72, 0.15, 15.1)),
  "data-green": pair(oklch(0.551, 0.153, 145.6), oklch(0.7, 0.15, 145.6)),
  "data-yellow": pair(oklch(0.785, 0.138, 101.1), oklch(0.8, 0.13, 101.1)),
  "data-cyan": pair(oklch(0.797, 0.106, 222.5), oklch(0.8, 0.1, 222.5)),
  "data-purple": pair(oklch(0.519, 0.169, 349.4), oklch(0.68, 0.15, 349.4)),
  "data-gray": pair(gray(0.792), gray(0.62)),
} as const satisfies Record<string, ThemePair>;

export const ALIAS_COLORS = {
  "card-foreground": "foreground",
  "popover-foreground": "foreground",
  secondary: "muted",
  "secondary-foreground": "foreground",
  accent: "muted",
  "accent-foreground": "foreground",
  sidebar: "card",
  "sidebar-foreground": "foreground",
  "sidebar-primary": "primary",
  "sidebar-primary-foreground": "primary-foreground",
  "sidebar-accent": "muted",
  "sidebar-accent-foreground": "foreground",
  "sidebar-border": "border",
  "sidebar-ring": "ring",
  "chart-paper": "card",
  "chart-ink": "foreground",
  "chart-muted": "muted-foreground",
  "chart-grid": "border",
  "chart-axis": "input",
  "chart-series": "data-blue",
  "chart-fit": "foreground",
  "chart-fault": "data-red",
  "chart-selection": "data-purple",
} as const satisfies Record<string, BaseColorName>;

export const DATA_COLORS = [
  "data-blue",
  "data-red",
  "data-green",
  "data-yellow",
  "data-cyan",
  "data-purple",
  "data-gray",
] as const satisfies readonly BaseColorName[];

export const CONTRAST_REQUIREMENTS: readonly ContrastRequirement[] = [
  ...pairs(["foreground"], ["background", "card", "popover", "muted", "selected"], 4.5),
  ...pairs(["muted-foreground"], ["background", "card", "popover", "muted", "selected"], 4.5),
  ...pairs(["link"], ["background", "card", "selected"], 4.5),
  ...pairs(["primary-foreground"], ["primary"], 4.5),
  ...pairs(["brand", "destructive", "warning", "success"], ["background", "card"], 4.5),
  ...pairs(["destructive"], ["destructive-subtle"], 4.5),
  ...pairs(["warning"], ["warning-subtle"], 4.5),
  ...pairs(["success"], ["success-subtle"], 4.5),
  ...pairs(["input", "ring"], ["background", "card"], 3),
  ...pairs(["data-blue", "data-red", "data-green", "data-purple"], ["card"], 3),
];

export const TYPEFACES = {
  sans: ['"IBM Plex Sans Variable"', '"IBM Plex Sans"', "ui-sans-serif", "system-ui", "sans-serif"],
  serif: ['"Source Serif 4 Variable"', '"Source Serif 4"', "Georgia", "serif"],
  mono: ['"IBM Plex Mono"', "ui-monospace", "monospace"],
} as const satisfies Record<string, readonly string[]>;

export const CORNERS = {
  sm: "3px",
  md: "4px",
  lg: "5px",
  xl: "8px",
} as const satisfies Record<string, string>;

export const MOTION = {
  duration: "120ms",
  easing: "cubic-bezier(0.2, 0, 0, 1)",
} as const satisfies Record<string, string>;

export function resolveColor(name: ColorName): ThemePair {
  return isAlias(name) ? BASE_COLORS[ALIAS_COLORS[name]] : BASE_COLORS[name];
}

export function isAlias(name: ColorName): name is AliasColorName {
  return Object.hasOwn(ALIAS_COLORS, name);
}

export type BaseColorName = keyof typeof BASE_COLORS;

export type AliasColorName = keyof typeof ALIAS_COLORS;

export type ColorName = BaseColorName | AliasColorName;

export interface ThemePair {
  light: Oklch;
  dark: Oklch;
}

export interface Oklch {
  l: number;
  c: number;
  h: number;
}

export interface ContrastRequirement {
  foreground: BaseColorName;
  background: BaseColorName;
  minimum: number;
}

function pairs(
  foregrounds: readonly BaseColorName[],
  backgrounds: readonly BaseColorName[],
  minimum: number,
): ContrastRequirement[] {
  return foregrounds.flatMap((foreground) => backgrounds.map((background) => ({ foreground, background, minimum })));
}

function pair(light: Oklch, dark: Oklch): ThemePair {
  return { light, dark };
}

function gray(l: number): Oklch {
  return oklch(l, 0, 0);
}

function blue(l: number, c: number): Oklch {
  return oklch(l, c, BLUE_HUE);
}

function oklch(l: number, c: number, h: number): Oklch {
  return { l, c, h };
}
