import type { IndexedFile } from "@neherlab/app-contracts";
import {
  Command,
  CommandCollection,
  CommandDialog,
  CommandEmpty,
  CommandGroup,
  CommandInput,
  CommandItem,
  CommandList,
} from "@neherlab/design/ui/command";
import { useCallback, useMemo, useState } from "react";

import { formatBytes } from "../format";
import { matchingExampleFiles } from "./exampleData";

export function ExampleFilesDialog({
  open,
  onOpenChange,
  title,
  files,
  choose,
}: {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  title: string;
  files: readonly IndexedFile[];
  choose: (file: IndexedFile) => void;
}) {
  return (
    <CommandDialog
      open={open}
      onOpenChange={onOpenChange}
      title={`Example files for ${title}`}
      description="Search the files of the example data"
      className="sm:max-w-2xl"
    >
      {open && <ExampleFileList files={files} choose={choose} />}
    </CommandDialog>
  );
}

function ExampleFileList({ files, choose }: { files: readonly IndexedFile[]; choose: (file: IndexedFile) => void }) {
  const [query, setQuery] = useState("");
  const groups = useMemo(() => [{ items: [...files] }], [files]);
  const filteredGroups = useMemo(() => [{ items: matchingExampleFiles(files, query) }], [files, query]);

  return (
    <Command
      items={groups}
      filteredItems={filteredGroups}
      value={query}
      onValueChange={setQuery}
      itemToStringValue={filePath}
    >
      <CommandInput placeholder="Search example files" aria-label="Search example files" />
      <CommandEmpty>Nothing matches.</CommandEmpty>
      <CommandList className="max-h-[min(60vh,32rem)]">
        {(group: { items: IndexedFile[] }) => (
          <CommandGroup key="files" items={group.items}>
            <CommandCollection>
              {(file: IndexedFile) => <ExampleFileEntry key={file.path} file={file} choose={choose} />}
            </CommandCollection>
          </CommandGroup>
        )}
      </CommandList>
    </Command>
  );
}

function ExampleFileEntry({ file, choose }: { file: IndexedFile; choose: (file: IndexedFile) => void }) {
  const onClick = useCallback(() => choose(file), [choose, file]);

  return (
    <CommandItem value={file} onClick={onClick}>
      <span className="min-w-0 flex-1 truncate font-mono text-xs">{file.path}</span>
      <span className="text-muted-foreground shrink-0 text-xs">{formatBytes(file.size)}</span>
    </CommandItem>
  );
}

function filePath(file: IndexedFile): string {
  return file.path;
}
