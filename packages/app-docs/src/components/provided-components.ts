import type { DesignSpecimen } from "@neherlab/design/specimen";
import type { Badge } from "@neherlab/design/ui/badge";
import type { Card, CardAction, CardDescription, CardHeader, CardTitle } from "@neherlab/design/ui/card";
import type { Kbd } from "@neherlab/design/ui/kbd";
import type { ReactElement, SVGProps } from "react";

import type { Disclosure } from "./Disclosure";
import type { FileTree } from "./FileTree";
import type { LinkButton } from "./LinkButton";
import type { CardGrid, LinkCard } from "./LinkCard";
import type { SyncedTabs } from "./SyncedTabs";

export interface ProvidedComponents {
  Badge: typeof Badge;
  Card: typeof Card;
  CardAction: typeof CardAction;
  CardDescription: typeof CardDescription;
  CardGrid: typeof CardGrid;
  CardHeader: typeof CardHeader;
  CardTitle: typeof CardTitle;
  DesignSpecimen: typeof DesignSpecimen;
  Disclosure: typeof Disclosure;
  FileTree: typeof FileTree;
  Kbd: typeof Kbd;
  LinkButton: typeof LinkButton;
  LinkCard: typeof LinkCard;
  SyncedTabs: typeof SyncedTabs;
}

type IconComponent = (props: SVGProps<SVGSVGElement>) => ReactElement;

// oxlint-disable-next-line custom/no-declare-global -- the MDX language server reads the components that pages use without an import only from the global type MDXProvidedComponents
declare global {
  type MDXProvidedComponents = ProvidedComponents & Record<`Icon${string}`, IconComponent>;
}
