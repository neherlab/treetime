export function ChartTooltip({ children }: { children: React.ReactNode }) {
  return (
    <div className="border-plate-border text-plate-ink rounded-md border bg-white px-2.5 py-1.5 text-xs shadow-sm">
      {children}
    </div>
  );
}
