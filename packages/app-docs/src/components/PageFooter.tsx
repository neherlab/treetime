import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { Card, CardDescription, CardTitle } from "@neherlab/app-ui/ui/card";
import ArrowLeft from "~icons/lucide/arrow-left";
import ArrowRight from "~icons/lucide/arrow-right";
import SquarePen from "~icons/lucide/square-pen";

export function PageFooter({
  previous,
  next,
  editUrl,
  locale,
}: {
  previous: FooterLink | undefined;
  next: FooterLink | undefined;
  editUrl: string | undefined;
  locale: Locale;
}) {
  return (
    <LocaleProvider locale={locale}>
      <footer className="not-prose mt-12 grid gap-6 border-t pt-6">
        <nav
          aria-label={`${m.docs_previous({}, { locale })} / ${m.docs_next({}, { locale })}`}
          className="grid gap-3 sm:grid-cols-2"
        >
          {previous !== undefined && (
            <NeighbourCard link={previous} label={m.docs_previous({}, { locale })} side="previous" />
          )}
          {next !== undefined && <NeighbourCard link={next} label={m.docs_next({}, { locale })} side="next" />}
        </nav>
        {editUrl !== undefined && (
          <a
            href={editUrl}
            className="text-muted-foreground hover:text-foreground inline-flex items-center gap-2 text-sm"
          >
            <SquarePen aria-hidden className="size-4" />
            {m.docs_edit_page({}, { locale })}
          </a>
        )}
      </footer>
    </LocaleProvider>
  );
}

export interface FooterLink {
  name: string;
  url: string;
}

function NeighbourCard({ link, label, side }: { link: FooterLink; label: string; side: "previous" | "next" }) {
  const Icon = side === "previous" ? ArrowLeft : ArrowRight;

  return (
    // oxlint-disable-next-line jsx-a11y/control-has-associated-label -- the card inside the link holds its text
    <a href={link.url} rel={side === "previous" ? "prev" : "next"} className={side === "next" ? "sm:col-start-2" : ""}>
      <Card size="sm" className="hover:bg-muted px-4 transition-colors">
        <CardDescription className="flex items-center gap-1.5">
          <Icon aria-hidden className="size-3.5" />
          {label}
        </CardDescription>
        <CardTitle>{link.name}</CardTitle>
      </Card>
    </a>
  );
}
