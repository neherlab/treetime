import { cn } from "cn";

import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { COLOR_GROUPS, colorFacts, type ColorRow, type WeakestContrast } from "./colors";

export function ColorSection() {
  return (
    <div className="flex flex-col gap-8">
      <p className="text-muted-foreground max-w-prose text-sm">
        Each token holds a light and a dark value. Grays have no hue, one blue marks what a person can act on, and lab
        red marks identity. The contrast column names the weakest pair that a test holds to its WCAG minimum in both
        themes: 4.5:1 for text, 3:1 for control outlines, focus, and chart marks.
      </p>
      {COLOR_GROUPS.map((group) => (
        <div key={group.title} className="flex flex-col gap-2">
          <h3 className="text-sm font-semibold">{group.title}</h3>
          <div className="overflow-x-auto rounded-lg border">
            <Table>
              <TableHeader>
                <TableRow>
                  <TableHead>Token</TableHead>
                  <TableHead>Use</TableHead>
                  <TableHead>Light</TableHead>
                  <TableHead>Dark</TableHead>
                  <TableHead>Weakest contrast</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {group.rows.map((row) => (
                  <ColorTableRow key={row.name} row={row} />
                ))}
              </TableBody>
            </Table>
          </div>
        </div>
      ))}
    </div>
  );
}

function ColorTableRow({ row }: { row: ColorRow }) {
  const facts = colorFacts(row.name);

  return (
    <TableRow>
      <TableCell className="font-mono text-xs">{row.name}</TableCell>
      <TableCell className="text-muted-foreground min-w-48 text-xs whitespace-normal">{row.use}</TableCell>
      <TableCell>
        <Swatch swatch={row.swatch} scheme="light" hex={facts.hex.light} />
      </TableCell>
      <TableCell>
        <Swatch swatch={row.swatch} scheme="dark" hex={facts.hex.dark} />
      </TableCell>
      <TableCell className="text-xs tabular-nums">
        {facts.weakest === undefined ? null : <ContrastText weakest={facts.weakest} />}
      </TableCell>
    </TableRow>
  );
}

function Swatch({ swatch, scheme, hex }: { swatch: string; scheme: "light" | "dark"; hex: string }) {
  return (
    <span className="flex items-center gap-2">
      <span
        aria-hidden
        className={cn(
          "inline-block size-6 shrink-0 rounded-md border",
          swatch,
          scheme === "light" ? "scheme-light" : "scheme-dark",
        )}
      />
      <span className="font-mono text-xs">{hex}</span>
    </span>
  );
}

function ContrastText({ weakest }: { weakest: WeakestContrast }) {
  const passes = Math.min(weakest.light, weakest.dark) >= weakest.minimum;

  return (
    <span className={passes ? undefined : "text-destructive"}>
      {weakest.light.toFixed(1)} and {weakest.dark.toFixed(1)}
      <span className="text-muted-foreground"> with {weakest.partner}</span>
    </span>
  );
}
