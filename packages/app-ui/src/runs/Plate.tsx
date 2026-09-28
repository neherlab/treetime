import type { ReactNode } from "react";

import { Card } from "../ui/card";

export function Plate({ title, caption, children }: { title: string; caption?: ReactNode; children: ReactNode }) {
  return (
    <Card size="sm" className="gap-0 py-0 shadow-none">
      <figure>
        <figcaption className="grid gap-0.5 border-b px-(--card-spacing) py-2.5">
          <span className="font-heading text-sm font-medium">{title}</span>
          {caption !== undefined && <span className="text-muted-foreground text-xs">{caption}</span>}
        </figcaption>
        <div className="p-2">{children}</div>
      </figure>
    </Card>
  );
}
