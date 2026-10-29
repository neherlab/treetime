import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

import IconsResolver from "unplugin-icons/resolver";

import type { ProvidedComponents } from "../components/provided-components";

const require = createRequire(import.meta.url);

const ICON_PREFIX = "Icon";

const resolveIcon = IconsResolver({ prefix: ICON_PREFIX, strict: true });

const PROVIDED_COMPONENTS = new Map<string, ProvidedComponent>(
  Object.entries({
    Badge: { module: "@neherlab/design/ui/badge" },
    Card: { module: "@neherlab/design/ui/card" },
    CardAction: { module: "@neherlab/design/ui/card" },
    CardDescription: { module: "@neherlab/design/ui/card" },
    CardGrid: { module: "../components/LinkCard.tsx" },
    CardHeader: { module: "@neherlab/design/ui/card" },
    CardTitle: { module: "@neherlab/design/ui/card" },
    DesignSpecimen: { module: "@neherlab/design/specimen", client: "visible" },
    Disclosure: { module: "../components/Disclosure.tsx", client: "visible" },
    FileTree: { module: "../components/FileTree.tsx" },
    Kbd: { module: "@neherlab/design/ui/kbd" },
    LinkButton: { module: "../components/LinkButton.tsx" },
    LinkCard: { module: "../components/LinkCard.tsx" },
    SyncedTabs: { module: "../components/SyncedTabs.tsx", client: "load" },
  } satisfies Record<keyof ProvidedComponents, ProvidedComponent>),
);

export function providedImport(name: string): ProvidedImport | undefined {
  const component = PROVIDED_COMPONENTS.get(name);

  if (component !== undefined) {
    return { source: modulePath(component.module), exported: name, client: component.client };
  }

  const icon = name.startsWith(ICON_PREFIX) ? resolveIcon(name) : undefined;

  return icon === undefined ? undefined : { source: icon, exported: "default", client: undefined };
}

export interface ProvidedImport {
  source: string;
  exported: string;
  client: Hydration | undefined;
}

interface ProvidedComponent {
  module: string;
  client?: Hydration;
}

type Hydration = "load" | "idle" | "visible";

function modulePath(module: string): string {
  return module.startsWith(".") ? fileURLToPath(new URL(module, import.meta.url)) : require.resolve(module);
}
