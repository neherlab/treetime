const FACES = [
  {
    name: "Source Serif 4",
    className: "font-serif",
    use: "Headings and page titles of sites and docs. Variable weight and optical size.",
  },
  {
    name: "IBM Plex Sans",
    className: "font-sans font-medium",
    use: "Interface and body text everywhere, and the text of charts. Tabular figures line up numbers in tables.",
  },
  {
    name: "IBM Plex Mono",
    className: "font-mono",
    use: "Code, commands, file names, node names, sequences, and wordmarks: text a person types or copies.",
  },
] as const;

const SCALE = [
  {
    px: 60,
    use: "Page title of a site",
    className: "font-serif text-6xl text-brand tracking-tight",
    sample: "Pathogen evolution",
  },
  { px: 48, use: "Article title", className: "font-serif text-5xl tracking-tight", sample: "Waves of immune escape" },
  {
    px: 36,
    use: "Page title of the docs; section",
    className: "font-serif text-4xl tracking-tight",
    sample: "Recent publications",
  },
  { px: 24, use: "Subsection", className: "font-serif text-2xl", sample: "Evolution of bacteria" },
  {
    px: 21,
    use: "Lead paragraph",
    className: "text-xl",
    sample: "TreeTime estimates the dates of ancestors from sampling dates.",
  },
  {
    px: 18,
    use: "Body of a site; app page title",
    className: "text-lg",
    sample: "The clock rate comes from a regression of divergence on date.",
  },
  {
    px: 16,
    use: "Body of the docs",
    className: "text-base",
    sample: "Set a clock rate to fix it instead of estimating it.",
  },
  {
    px: 14,
    use: "Body of an app",
    className: "text-sm",
    sample: "Clock rate 3.54 ± 0.21 × 10⁻³ substitutions per site per year",
  },
  {
    px: 12,
    use: "Captions, table headers, hints",
    className: "text-xs text-muted-foreground",
    sample: "90% interval from the marginal time posterior",
  },
] as const;

export function TypeSection() {
  return (
    <div className="flex flex-col gap-8">
      <div className="grid gap-4 @3xl:grid-cols-3">
        {FACES.map((face) => (
          <div key={face.name} className="bg-card flex flex-col gap-2 rounded-lg border p-5">
            <span aria-hidden className={`${face.className} text-5xl leading-none`}>
              Ag
            </span>
            <span className="text-sm font-semibold">{face.name}</span>
            <span className="text-muted-foreground text-sm">{face.use}</span>
          </div>
        ))}
      </div>
      <div className="divide-y border-y">
        {SCALE.map((step) => (
          <div
            key={step.px}
            className="grid grid-cols-[3rem_minmax(0,1fr)] items-baseline gap-x-4 py-3 @3xl:grid-cols-[3rem_14rem_minmax(0,1fr)]"
          >
            <span className="text-muted-foreground font-mono text-xs">{step.px}</span>
            <span className="text-muted-foreground hidden text-xs @3xl:block">{step.use}</span>
            <span className={`${step.className} truncate`}>{step.sample}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
