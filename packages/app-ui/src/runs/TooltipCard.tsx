export function TooltipCard({ children }: { children: React.ReactNode }) {
  return (
    <div className="border-border/50 bg-background grid min-w-32 items-start gap-1 rounded-lg border px-2.5 py-1.5 text-xs shadow-xl">
      {children}
    </div>
  );
}
