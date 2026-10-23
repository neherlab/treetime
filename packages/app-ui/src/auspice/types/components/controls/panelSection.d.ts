import type { ComponentType, ReactElement } from "react";

export declare const PanelSection: ComponentType<{
  panel: "tree" | "map" | "entropy";
  title: string;
  tooltip?: ReactElement;
  options?: ReactElement;
}>;
