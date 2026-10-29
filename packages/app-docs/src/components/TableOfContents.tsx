import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import { cn } from "cn";
import { AnchorProvider, useActiveAnchor } from "fumadocs-core/toc";
import { useMemo } from "react";

export function TableOfContents({ headings, locale }: { headings: TocHeading[]; locale: Locale }) {
  const toc = useMemo(
    () => headings.map((heading) => ({ title: heading.text, url: `#${heading.slug}`, depth: heading.depth })),
    [headings],
  );

  return (
    <LocaleProvider locale={locale}>
      <AnchorProvider toc={toc}>
        <nav aria-labelledby="toc-title" className="text-sm">
          <p id="toc-title" className="mb-2 font-semibold">
            {m.docs_toc({}, { locale })}
          </p>
          <TocList headings={headings} />
        </nav>
      </AnchorProvider>
    </LocaleProvider>
  );
}

export interface TocHeading {
  depth: number;
  slug: string;
  text: string;
}

function TocList({ headings }: { headings: TocHeading[] }) {
  const active = useActiveAnchor();

  return (
    <ul className="grid gap-1.5 border-l">
      {headings.map((heading) => (
        <li key={heading.slug}>
          <a
            href={`#${heading.slug}`}
            aria-current={heading.slug === active ? "location" : undefined}
            className={cn(
              "text-muted-foreground hover:text-foreground aria-[current=location]:border-primary aria-[current=location]:text-foreground -ml-px block border-l border-transparent py-0.5 pl-3",
              heading.depth === 3 && "pl-6",
            )}
          >
            {heading.text}
          </a>
        </li>
      ))}
    </ul>
  );
}
