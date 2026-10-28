import { SearchTrigger } from "@neherlab/app-ui/chrome/SearchTrigger";
import { PALETTE_HOTKEY } from "@neherlab/app-ui/hotkeys";
import { m, useLocale } from "@neherlab/app-ui/i18n";
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
import { useHotkey } from "@tanstack/react-hotkeys";
import { useCallback, useEffect, useMemo, useState } from "react";

import { excerptParts, loadPagefind, type SearchPages, type SearchResult } from "../lib/pagefind";

export function DocsSearch({ loadSearch = loadPagefind }: { loadSearch?: () => Promise<SearchPages> }) {
  const locale = useLocale();
  const [open, setOpen] = useState(false);
  const show = useCallback(() => setOpen(true), []);

  useHotkey(PALETTE_HOTKEY, show);
  useHotkey("/", show);

  return (
    <>
      <SearchTrigger label={m.docs_search({}, { locale })} hotkey={PALETTE_HOTKEY} onClick={show} />
      <CommandDialog
        open={open}
        onOpenChange={setOpen}
        title={m.docs_search({}, { locale })}
        description={m.docs_search_placeholder({}, { locale })}
        className="sm:max-w-2xl"
      >
        {open && <SearchBody loadSearch={loadSearch} />}
      </CommandDialog>
    </>
  );
}

function SearchBody({ loadSearch }: { loadSearch: () => Promise<SearchPages> }) {
  const locale = useLocale();
  const search = useSearch(loadSearch);
  const [query, setQuery] = useState("");
  const results = useResults(search, query);
  const groups = useMemo(() => [{ items: results }], [results]);

  return (
    <Command
      items={groups}
      filteredItems={groups}
      value={query}
      onValueChange={setQuery}
      itemToStringValue={resultTitle}
    >
      <CommandInput
        placeholder={m.docs_search_placeholder({}, { locale })}
        aria-label={m.docs_search({}, { locale })}
      />
      {search.kind === "failed" ? (
        <p className="text-muted-foreground px-4 py-6 text-center text-sm">
          {m.docs_search_unavailable({}, { locale })}
        </p>
      ) : (
        <>
          {query.trim() !== "" && search.kind === "ready" && (
            <CommandEmpty>{m.docs_search_empty({}, { locale })}</CommandEmpty>
          )}
          {search.kind === "loading" && (
            <p className="text-muted-foreground px-4 py-6 text-center text-sm">
              {m.docs_search_loading({}, { locale })}
            </p>
          )}
          <CommandList className="max-h-[min(60vh,32rem)]" aria-label={m.docs_search_results({}, { locale })}>
            {(group: { items: SearchResult[] }) => (
              <CommandGroup key="results" items={group.items}>
                <CommandCollection>
                  {(result: SearchResult) => <ResultEntry key={result.url} result={result} />}
                </CommandCollection>
              </CommandGroup>
            )}
          </CommandList>
        </>
      )}
    </Command>
  );
}

function ResultEntry({ result }: { result: SearchResult }) {
  const open = useCallback(() => window.location.assign(result.url), [result.url]);

  return (
    <CommandItem value={result} onClick={open} className="items-start">
      <a href={result.url} className="grid min-w-0 gap-1" tabIndex={-1}>
        <span className="font-medium">{resultTitle(result)}</span>
        <Excerpt text={result.excerpt} />
        {result.sub_results
          .filter((sub) => sub.url !== result.url)
          .slice(0, 3)
          .map((sub) => (
            <span key={sub.url} className="border-l pl-2 text-xs">
              <span className="text-foreground">{sub.title}</span> <Excerpt text={sub.excerpt} />
            </span>
          ))}
      </a>
    </CommandItem>
  );
}

function Excerpt({ text }: { text: string }) {
  return (
    <span className="text-muted-foreground line-clamp-2 text-xs">
      {excerptParts(text).map((part, index) =>
        part.marked ? (
          // oxlint-disable-next-line react/no-array-index-key -- excerpt parts have no identity beyond their position
          <mark key={index} className="bg-accent text-accent-foreground rounded-sm">
            {part.text}
          </mark>
        ) : (
          part.text
        ),
      )}
    </span>
  );
}

type SearchState = { kind: "loading" } | { kind: "ready"; search: SearchPages } | { kind: "failed" };

function useSearch(loadSearch: () => Promise<SearchPages>): SearchState {
  const [state, setState] = useState<SearchState>({ kind: "loading" });

  useEffect(() => {
    let active = true;

    loadSearch().then(
      (search) => active && setState({ kind: "ready", search }),
      () => active && setState({ kind: "failed" }),
    );

    return () => {
      active = false;
    };
  }, [loadSearch]);

  return state;
}

function useResults(search: SearchState, query: string): SearchResult[] {
  const [results, setResults] = useState<SearchResult[]>([]);

  useEffect(() => {
    if (search.kind !== "ready" || query.trim() === "") {
      return undefined;
    }

    let active = true;

    void search.search(query).then((found) => active && setResults(found));

    return () => {
      active = false;
    };
  }, [query, search]);

  return search.kind === "ready" && query.trim() !== "" ? results : [];
}

function resultTitle(result: SearchResult): string {
  return result.meta.title ?? result.url;
}
