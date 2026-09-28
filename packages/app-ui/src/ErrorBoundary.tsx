import { errorMessage } from "@neherlab/app-contracts";
import { ErrorBoundary as ReactErrorBoundary, type FallbackProps } from "react-error-boundary";

import { Button } from "./ui/button";
import { Empty, EmptyContent, EmptyDescription, EmptyHeader, EmptyTitle } from "./ui/empty";

function ErrorFallback({ error, resetErrorBoundary }: FallbackProps) {
  return (
    <Empty role="alert" className="min-h-svh">
      <EmptyHeader>
        <EmptyTitle>Something went wrong</EmptyTitle>
        <EmptyDescription>The app stopped because of an error it could not handle.</EmptyDescription>
      </EmptyHeader>
      <EmptyContent className="max-w-2xl">
        <pre className="bg-muted text-muted-foreground max-h-64 w-full overflow-auto rounded-md p-4 text-left font-mono text-xs">
          {errorMessage(error)}
        </pre>
        <Button variant="outline" onClick={resetErrorBoundary}>
          Try again
        </Button>
      </EmptyContent>
    </Empty>
  );
}

export function ErrorBoundary({ children }: { children: React.ReactNode }) {
  return <ReactErrorBoundary FallbackComponent={ErrorFallback}>{children}</ReactErrorBoundary>;
}
