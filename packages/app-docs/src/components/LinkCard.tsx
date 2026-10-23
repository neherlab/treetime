import { Card, CardDescription, CardHeader, CardTitle } from "@neherlab/app-ui/ui/card";
import ArrowRight from "~icons/lucide/arrow-right";

export function LinkCard({ title, description, href }: { title: string; description?: string; href: string }) {
  return (
    // oxlint-disable-next-line jsx-a11y/control-has-associated-label -- the card inside the link holds its text
    <a href={href} className="not-prose group/link block">
      <Card size="sm" className="group-hover/link:bg-muted h-full transition-colors">
        <CardHeader>
          <CardTitle className="flex items-center justify-between gap-2">
            {title}
            <ArrowRight aria-hidden className="text-muted-foreground size-4" />
          </CardTitle>
          {description !== undefined && <CardDescription>{description}</CardDescription>}
        </CardHeader>
      </Card>
    </a>
  );
}

export function CardGrid({ children }: { children?: React.ReactNode }) {
  return <div className="not-prose my-6 grid gap-3 sm:grid-cols-2">{children}</div>;
}
