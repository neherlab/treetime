export function ChartTooltip({ children }: { children: React.ReactNode }) {
  return (
    <div className="border-line-strong text-ink bg-surface-1 rounded-md border px-2.5 py-1.5 text-xs shadow-sm">
      {children}
    </div>
  );
}
