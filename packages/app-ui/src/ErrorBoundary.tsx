import { ErrorBoundary as ReactErrorBoundary, type FallbackProps } from "react-error-boundary";

import { Button } from "./ui";

function ErrorFallback({ error, resetErrorBoundary }: FallbackProps) {
  const message = error instanceof Error ? error.message : String(error);
  return (
    <div
      role="alert"
      className="bg-surface-0 text-ink flex h-screen flex-col items-center justify-center gap-4 p-8"
    >
      <h1 className="text-lg font-semibold">Something went wrong</h1>
      <pre className="bg-surface-1 text-ink-muted max-h-64 max-w-2xl overflow-auto rounded-md p-4 text-sm">
        {message}
      </pre>
      <Button variant="outline" onClick={resetErrorBoundary}>
        Try again
      </Button>
    </div>
  );
}

export function ErrorBoundary({ children }: { children: React.ReactNode }) {
  return <ReactErrorBoundary FallbackComponent={ErrorFallback}>{children}</ReactErrorBoundary>;
}
