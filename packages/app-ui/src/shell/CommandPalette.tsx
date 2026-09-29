import { errorMessage, type ExampleConfig } from "@neherlab/app-contracts";
import { datasets, runsList } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useTheme } from "next-themes";
import { useCallback, useMemo, useRef, useState } from "react";

import { settingFieldId } from "../analysis/fieldIds";
import { useConfigLoader } from "../analysis/useConfigLoader";
import { useApi } from "../api/hooks";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { useDraftStore } from "../store/draft";
import { useShellStore } from "../store/shell";
import {
  Command,
  CommandCollection,
  CommandDialog,
  CommandEmpty,
  CommandGroup,
  CommandGroupLabel,
  CommandInput,
  CommandItem,
  CommandList,
} from "../ui/command";
import { useToastManager } from "../ui/toast";
import { matchingPaletteItems, type PaletteGroup, type PaletteItem, paletteGroups, paletteItem } from "./paletteItems";
import { changedFlags, listedRuns } from "./runList";
import { nextTheme } from "./SiteHeader";

export function CommandPalette() {
  const open = useShellStore((state) => state.paletteOpen);
  const setOpen = useShellStore((state) => state.setPaletteOpen);
  const focusTarget = useRef<string | undefined>(undefined);

  const onOpenChange = useCallback(
    (next: boolean) => {
      if (next) {
        focusTarget.current = undefined;
      }

      setOpen(next);
    },
    [setOpen],
  );

  const finalFocus = useCallback(
    () => (focusTarget.current === undefined ? true : (document.getElementById(focusTarget.current) ?? true)),
    [],
  );

  const choose = useCallback(
    async (item: PaletteItem) => {
      focusTarget.current = item.focusId;

      if (item.focusId === undefined) {
        setOpen(false);
        await item.run();
      } else {
        await item.run();
        setOpen(false);
      }
    },
    [setOpen],
  );

  return (
    <CommandDialog
      open={open}
      onOpenChange={onOpenChange}
      title="Search"
      description="Search runs, settings, examples and actions"
      className="sm:max-w-2xl"
      finalFocus={finalFocus}
    >
      {open && <PaletteBody choose={choose} />}
    </CommandDialog>
  );
}

function PaletteBody({ choose }: { choose: (item: PaletteItem) => Promise<void> }) {
  const items = usePaletteItems();
  const [query, setQuery] = useState("");
  const groups = useMemo(() => paletteGroups(items), [items]);
  const filteredGroups = useMemo(() => paletteGroups(matchingPaletteItems(items, query)), [items, query]);

  return (
    <Command
      items={groups}
      filteredItems={filteredGroups}
      value={query}
      onValueChange={setQuery}
      itemToStringValue={itemTitle}
    >
      <CommandInput placeholder="Search runs, settings, examples" aria-label="Search" />
      <CommandEmpty>Nothing matches.</CommandEmpty>
      <CommandList className="max-h-[min(60vh,32rem)]">
        {(group: PaletteGroup) => (
          <CommandGroup key={group.kind} items={group.items}>
            <CommandGroupLabel>{`${group.kind}s`}</CommandGroupLabel>
            <CommandCollection>
              {(item: PaletteItem) => <PaletteEntry key={item.id} item={item} choose={choose} />}
            </CommandCollection>
          </CommandGroup>
        )}
      </CommandList>
    </Command>
  );
}

function PaletteEntry({ item, choose }: { item: PaletteItem; choose: (item: PaletteItem) => Promise<void> }) {
  const onClick = useCallback(() => void choose(item), [choose, item]);

  return (
    <CommandItem value={item} onClick={onClick} className="items-start">
      <span className="grid min-w-0 gap-0.5">
        <span className="truncate">{item.title}</span>
        {item.description !== "" && (
          <span className="text-muted-foreground line-clamp-1 text-xs">{item.description}</span>
        )}
      </span>
    </CommandItem>
  );
}

function itemTitle(item: PaletteItem): string {
  return item.title;
}

function usePaletteItems(): PaletteItem[] {
  const navigate = useNavigate();
  const { data: runList } = useApi((context) => runsList(context));
  const { data: catalog } = useApi((context) => datasets(context), { staleTime: Infinity });
  const command = useDraftStore((state) => state.command);
  const compareIds = useShellStore((state) => state.compareIds);
  const { theme, setTheme } = useTheme();
  const loadConfig = useConfigLoader();
  const toasts = useToastManager();

  const loadExample = useCallback(
    async (example: ExampleConfig) => {
      try {
        const result = await loadConfig(example.content, example.command, false);

        if (!result.loaded) {
          toasts.add({ title: `${example.path} cannot be loaded`, description: result.messages.join("; ") });
        }
      } catch (error: unknown) {
        toasts.add({ title: `${example.path} cannot be loaded`, description: errorMessage(error) });
      }
    },
    [loadConfig, toasts],
  );

  return useMemo(() => {
    const items: PaletteItem[] = [
      paletteItem("Action", "action-new", "New analysis", "", () => navigate({ to: "/new" })),
      paletteItem("Action", "action-theme", "Change theme", "System, light or dark", () => setTheme(nextTheme(theme))),
    ];

    const [first, second] = compareIds;

    if (first !== undefined && second !== undefined) {
      items.push(
        paletteItem("Action", "action-compare", "Compare selected runs", "", () =>
          navigate({ to: "/compare/$a/$b", params: { a: first, b: second } }),
        ),
      );
    }

    for (const run of listedRuns(runList?.runs ?? [], "", null)) {
      items.push(
        paletteItem(
          "Run",
          `run-${run.id}`,
          run.title,
          [COMMAND_INFO[run.command].label, ...changedFlags(run)].join("  "),
          () => navigate({ to: "/runs/$id/results", params: { id: run.id } }),
        ),
      );
    }

    for (const spec of COMMAND_SETTINGS[command].specs) {
      if (spec.role !== "output") {
        items.push(
          paletteItem(
            "Setting",
            `setting-${spec.key}`,
            `${spec.label}  ${spec.flag}`,
            spec.help,
            () => {
              useDraftStore.getState().update({ view: "all", search: spec.key, changedOnly: false });

              return navigate({ to: "/new" });
            },
            settingFieldId(spec.key),
          ),
        );
      }
    }

    for (const example of catalog?.examples ?? []) {
      items.push(
        paletteItem("Example", `example-${example.path}`, example.title, example.path, () => loadExample(example)),
      );
    }

    return items;
  }, [catalog, command, compareIds, loadExample, navigate, runList, setTheme, theme]);
}
