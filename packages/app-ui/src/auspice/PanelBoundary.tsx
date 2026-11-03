import { errorMessage } from "@neherlab/app-contracts";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { useCallback, type ReactNode } from "react";
import { ErrorBoundary, type FallbackProps } from "react-error-boundary";

import type { ViewPanel } from "./layout";

const FAILURE_TITLES: Readonly<Record<ViewPanel, string>> = {
  tree: "Auspice cannot draw the tree",
  map: "Auspice cannot draw the map",
  entropy: "Auspice cannot draw the entropy chart",
};

export function PanelBoundary({
  panel,
  onFailure,
  children,
}: {
  panel: ViewPanel;
  onFailure: PanelFailureHandler;
  children: ReactNode;
}) {
  const failed = useCallback(() => onFailure(panel), [onFailure, panel]);

  const fallback = useCallback(
    ({ error }: FallbackProps) => <PanelFailure title={FAILURE_TITLES[panel]} error={error} />,
    [panel],
  );

  return (
    <ErrorBoundary fallbackRender={fallback} onError={failed}>
      {children}
    </ErrorBoundary>
  );
}

export type PanelFailureHandler = (panel: ViewPanel) => void;

function PanelFailure({ title, error }: { title: string; error: unknown }) {
  return (
    <Alert variant="destructive">
      <AlertTitle>{title}</AlertTitle>
      <AlertDescription>{errorMessage(error)}</AlertDescription>
    </Alert>
  );
}
