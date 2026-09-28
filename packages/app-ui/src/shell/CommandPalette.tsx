import { errorMessage, type ExampleConfig } from "@neherlab/app-contracts";
import { datasets, runsList } from "@neherlab/app-contracts/client";
import { useNavigate } from "@tanstack/react-router";
import { useTheme } from "next-themes";
import { useCallback, useMemo } from "react";

import { settingFieldId } from "../analysis/fieldIds";
import { useConfigLoader } from "../analysis/useConfigLoader";
import { useApi } from "../api/hooks";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { useDraftStore } from "../store/draft";
import { useShellStore } from "../store/shell";
import {
  Command,
  CommandDialog,
  CommandEmpty,
  CommandGroup,
  CommandInput,
  CommandItem,
  CommandList,
} from "../ui/command";
import { useToastManager } from "../ui/toast";
import { changedFlags, listedRuns } from "./runList";
import { nextTheme } from "./SiteHeader";

const PALETTE_KINDS = ["Action", "Run", "Setting", "Example"] as const;

type PaletteKind = (typeof PALETTE_KINDS)[number];

const FOCUS_DELAY_MS = 50;

export function CommandPalette() {
  const open = useShellStore((state) => state.paletteOpen);
  const setOpen = useShellStore((state) => state.setPaletteOpen);
  const close = useCallback(() => setOpen(false), [setOpen]);

  return (
    <CommandDialog
      open={open}
      onOpenChange={setOpen}
      title="Search"
      description="Search runs, settings, examples and actions"
      className="sm:max-w-2xl"
    >
      {open && <PaletteBody close={close} />}
    </CommandDialog>
  );
}

function PaletteBody({ close }: { close: () => void }) {
  const items = usePaletteItems();

  return (
    <Command>
      <CommandInput placeholder="Search runs, settings, examples" />
      <CommandList className="max-h-[min(60vh,32rem)]">
        <CommandEmpty>Nothing matches.</CommandEmpty>
        {PALETTE_KINDS.map((kind) => (
          <PaletteGroup key={kind} kind={kind} items={items} close={close} />
        ))}
      </CommandList>
    </Command>
  );
}

function PaletteGroup({ kind, items, close }: { kind: PaletteKind; items: readonly PaletteItem[]; close: () => void }) {
  const members = items.filter((item) => item.kind === kind);

  if (members.length === 0) {
    return null;
  }

  return (
    <CommandGroup heading={`${kind}s`}>
      {members.map((item) => (
        <PaletteEntry key={item.id} item={item} close={close} />
      ))}
    </CommandGroup>
  );
}

function PaletteEntry({ item, close }: { item: PaletteItem; close: () => void }) {
  const onSelect = useCallback(() => {
    close();
    item.run();
  }, [close, item]);

  return (
    <CommandItem value={item.id} keywords={item.keywords} onSelect={onSelect} className="items-start">
      <span className="grid min-w-0 gap-0.5">
        <span className="truncate">{item.title}</span>
        {item.description !== "" && (
          <span className="text-muted-foreground line-clamp-1 text-xs">{item.description}</span>
        )}
      </span>
    </CommandItem>
  );
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
      paletteItem("Action", "action-new", "New analysis", "", () => void navigate({ to: "/new" })),
      paletteItem("Action", "action-theme", "Change theme", "System, light or dark", () => setTheme(nextTheme(theme))),
    ];

    const [first, second] = compareIds;

    if (first !== undefined && second !== undefined) {
      items.push(
        paletteItem(
          "Action",
          "action-compare",
          "Compare selected runs",
          "",
          () => void navigate({ to: "/compare/$a/$b", params: { a: first, b: second } }),
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
          () => void navigate({ to: "/runs/$id/results", params: { id: run.id } }),
        ),
      );
    }

    for (const spec of COMMAND_SETTINGS[command].specs) {
      if (spec.role !== "output") {
        items.push(
          paletteItem("Setting", `setting-${spec.key}`, `${spec.label}  ${spec.flag}`, spec.help, () => {
            useDraftStore.getState().update({ view: "all", search: spec.key, changedOnly: false });
            void navigate({ to: "/new" }).then(() => focusLater(settingFieldId(spec.key)));
          }),
        );
      }
    }

    for (const example of catalog?.examples ?? []) {
      items.push(
        paletteItem("Example", `example-${example.path}`, example.title, example.path, () => void loadExample(example)),
      );
    }

    return items;
  }, [catalog, command, compareIds, loadExample, navigate, runList, setTheme, theme]);
}

function paletteItem(kind: PaletteKind, id: string, title: string, description: string, run: () => void): PaletteItem {
  return { id, kind, title, description, keywords: [kind, title, description], run };
}

function focusLater(id: string) {
  globalThis.setTimeout(() => document.getElementById(id)?.focus(), FOCUS_DELAY_MS);
}

interface PaletteItem {
  id: string;
  kind: PaletteKind;
  title: string;
  description: string;
  keywords: string[];
  run: () => void;
}
