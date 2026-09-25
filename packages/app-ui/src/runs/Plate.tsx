import type { ReactNode } from "react";

export function Plate({ title, caption, children }: { title: string; caption?: ReactNode; children: ReactNode }) {
  return (
    <figure className="border-line bg-surface-1 m-0 overflow-hidden rounded-lg border">
      <figcaption className="border-line flex flex-wrap items-baseline gap-x-2.5 gap-y-0.5 border-b px-3.5 py-2.5">
        <span className="font-bold">{title}</span>
        {caption !== undefined && <span className="text-ink-faint text-xs">{caption}</span>}
      </figcaption>
      <div className="plate px-2 py-2">{children}</div>
    </figure>
  );
}
