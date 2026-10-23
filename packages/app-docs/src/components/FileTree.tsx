import File from "~icons/lucide/file";
import Folder from "~icons/lucide/folder";

export function FileTree({ entries }: { entries: FileTreeEntry[] }) {
  return (
    <div className="not-prose bg-card my-6 rounded-lg border p-3 font-mono text-sm">
      <TreeLevel entries={entries} />
    </div>
  );
}

export interface FileTreeEntry {
  name: string;
  note?: string;
  children?: FileTreeEntry[];
}

function TreeLevel({ entries }: { entries: FileTreeEntry[] }) {
  return (
    <ul className="grid gap-1 [&_ul]:mt-1 [&_ul]:ml-4 [&_ul]:border-l [&_ul]:pl-3">
      {entries.map((entry) => (
        <li key={entry.name}>
          <span className="inline-flex items-center gap-2">
            {entry.children === undefined ? (
              <File aria-hidden className="text-muted-foreground size-4" />
            ) : (
              <Folder aria-hidden className="text-primary size-4" />
            )}
            {entry.name}
            {entry.note !== undefined && <span className="text-muted-foreground font-sans text-xs">{entry.note}</span>}
          </span>
          {entry.children !== undefined && <TreeLevel entries={entry.children} />}
        </li>
      ))}
    </ul>
  );
}
