import { Wordmark } from "../ui/wordmark";

const WORDMARKS = [
  { name: "neherlab", place: "biozentrum", use: "Lab site, talks, course pages" },
  { name: "treetime", place: undefined, use: "TreeTime app, CLI, and docs" },
  { name: "treeknit", place: undefined, use: "TreeKnit app" },
  { name: "pangraph", place: undefined, use: "Pangraph docs" },
] as const;

export function IdentitySection() {
  return (
    <div className="grid gap-8 @3xl:grid-cols-[3fr_2fr]">
      <ul className="divide-y border-y">
        {WORDMARKS.map(({ name, place, use }) => (
          <li key={name} className="flex items-baseline justify-between gap-4 py-4">
            {place === undefined ? (
              <Wordmark name={name} className="text-2xl" />
            ) : (
              <Wordmark name={name} place={place} className="text-2xl" />
            )}
            <span className="text-muted-foreground text-right text-xs">{use}</span>
          </li>
        ))}
      </ul>
      <div className="flex flex-col gap-4 text-sm">
        <p>
          The wordmark of the lab site reads like a shell prompt: a name, an <code>@</code>, a place. A tool shows its
          name in the same monospace face, which is also the command that runs it.
        </p>
        <p className="text-muted-foreground">
          Logos follow the lab mark: lines in gray, at most one stroke in coral and one tip in blue.
        </p>
        <div className="flex gap-2" aria-hidden>
          <span className="bg-mark-gray h-2 w-10 rounded-sm" />
          <span className="bg-mark-coral h-2 w-10 rounded-sm" />
          <span className="bg-mark-blue h-2 w-10 rounded-sm" />
        </div>
      </div>
    </div>
  );
}
