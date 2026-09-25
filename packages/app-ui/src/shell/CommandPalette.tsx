import { useNavigate } from "@tanstack/react-router";
import { useTheme } from "next-themes";
import { useCallback, useMemo, useRef, useState } from "react";

import { settingFieldId } from "../analysis/fieldIds";
import { useConfigLoader } from "../analysis/useConfigLoader";
import { useDatasetCatalog, useRunList } from "../queries";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO, isAppCommand } from "../settings/commands";
import { settingLabel } from "../settings/labels";
import { useDraftStore } from "../store/draft";
import { useShellStore } from "../store/shell";
import { Dialog, Toast, cn } from "../ui";
import { changedFlags, listedRuns } from "./runList";
import { nextTheme } from "./TopBar";

interface PaletteItem {
  id: string;
  kind: "Action" | "Run" | "Setting" | "Example";
  title: string;
  description: string;
  run: () => void;
}

const ITEMS_SHOWN = 60;

const DESCRIPTION_LENGTH = 120;

const FOCUS_DELAY_MS = 50;

export function CommandPalette() {
  const open = useShellStore((state) => state.paletteOpen);
  const setOpen = useShellStore((state) => state.setPaletteOpen);
  const input = useRef<HTMLInputElement>(null);
  const close = useCallback(() => setOpen(false), [setOpen]);

  return (
    <Dialog.Root open={open} onOpenChange={setOpen}>
      <Dialog.Popup initialFocus={input} className="top-[12vh] max-w-2xl translate-y-0 p-0" aria-label="Search">
        {open && <PaletteBody input={input} close={close} />}
      </Dialog.Popup>
    </Dialog.Root>
  );
}

function PaletteBody({ input, close }: { input: React.RefObject<HTMLInputElement>; close: () => void }) {
  const items = usePaletteItems();
  const [query, setQuery] = useState("");
  const [active, setActive] = useState(0);

  const shown = useMemo(() => {
    const words = query
      .toLowerCase()
      .split(/\s+/u)
      .filter((word) => word !== "");

    return items
      .filter((item) => {
        const haystack = `${item.kind} ${item.title} ${item.description}`.toLowerCase();

        return words.every((word) => haystack.includes(word));
      })
      .slice(0, ITEMS_SHOWN);
  }, [items, query]);

  const choose = useCallback(
    (item: PaletteItem | undefined) => {
      if (item !== undefined) {
        close();
        item.run();
      }
    },
    [close],
  );

  const onKeyDown = useCallback(
    (event: React.KeyboardEvent) => {
      if (event.key === "ArrowDown") {
        event.preventDefault();
        setActive(Math.min(shown.length - 1, active + 1));
      } else if (event.key === "ArrowUp") {
        event.preventDefault();
        setActive(Math.max(0, active - 1));
      } else if (event.key === "Enter") {
        event.preventDefault();
        choose(shown[active]);
      }
    },
    [active, choose, shown],
  );

  const onQuery = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setQuery(event.target.value);
    setActive(0);
  }, []);

  return (
    <div className="grid">
      <input
        ref={input}
        type="text"
        value={query}
        onChange={onQuery}
        onKeyDown={onKeyDown}
        placeholder="Search runs, settings, examples"
        aria-label="Search"
        className="border-line w-full border-b bg-transparent px-4 py-3.5 text-base outline-none"
      />
      <ul aria-label="Results" className="max-h-[50vh] overflow-auto p-1.5">
        {shown.map((item, index) => (
          <PaletteEntry
            key={item.id}
            item={item}
            active={index === active}
            index={index}
            choose={choose}
            hover={setActive}
          />
        ))}
        {shown.length === 0 && <li className="text-ink-muted p-4 text-center">Nothing matches.</li>}
      </ul>
    </div>
  );
}

function PaletteEntry({
  item,
  active,
  index,
  choose,
  hover,
}: {
  item: PaletteItem;
  active: boolean;
  index: number;
  choose: (item: PaletteItem) => void;
  hover: (index: number) => void;
}) {
  const onClick = useCallback(() => choose(item), [choose, item]);
  const onMouseMove = useCallback(() => hover(index), [hover, index]);

  const description =
    item.description.length > DESCRIPTION_LENGTH
      ? `${item.description.slice(0, DESCRIPTION_LENGTH - 3)}...`
      : item.description;

  return (
    <li>
      <button
        type="button"
        tabIndex={-1}
        aria-current={active}
        onClick={onClick}
        onMouseMove={onMouseMove}
        className={cn(
          "grid w-full cursor-pointer grid-cols-[4.5rem_1fr] gap-x-2.5 rounded-md px-2.5 py-1.5 text-left",
          active && "bg-accent-subtle",
        )}
      >
        <span className="text-ink-faint text-xs">{item.kind}</span>
        <span>{item.title}</span>
        {description !== "" && <span className="text-ink-faint col-start-2 text-xs">{description}</span>}
      </button>
    </li>
  );
}

function usePaletteItems(): PaletteItem[] {
  const navigate = useNavigate();
  const { data: runList } = useRunList();
  const { data: catalog } = useDatasetCatalog();
  const command = useDraftStore((state) => state.command);
  const compareIds = useShellStore((state) => state.compareIds);
  const { theme, setTheme } = useTheme();
  const loadConfig = useConfigLoader();
  const toasts = Toast.useToastManager();

  return useMemo(() => {
    const items: PaletteItem[] = [
      {
        id: "action-new",
        kind: "Action",
        title: "New analysis",
        description: "",
        run: () => void navigate({ to: "/new" }),
      },
      {
        id: "action-theme",
        kind: "Action",
        title: "Change theme",
        description: "System, light or dark",
        run: () => setTheme(nextTheme(theme)),
      },
    ];

    const [first, second] = compareIds;

    if (first !== undefined && second !== undefined) {
      items.push({
        id: "action-compare",
        kind: "Action",
        title: "Compare selected runs",
        description: "",
        run: () => void navigate({ to: "/compare/$a/$b", params: { a: first, b: second } }),
      });
    }

    for (const run of listedRuns(runList?.runs ?? [], "", null)) {
      items.push({
        id: `run-${run.id}`,
        kind: "Run",
        title: run.title,
        description: [COMMAND_INFO[run.command].label, ...changedFlags(run)].join("  "),
        run: () => void navigate({ to: "/runs/$id/results", params: { id: run.id } }),
      });
    }

    for (const spec of COMMAND_SETTINGS[command].specs) {
      if (spec.pathRole === "output") {
        continue;
      }

      items.push({
        id: `setting-${spec.key}`,
        kind: "Setting",
        title: `${settingLabel(spec.key)}  ${spec.flag}`,
        description: spec.help,
        run: () => {
          useDraftStore.getState().update({ view: "all", search: spec.key, changedOnly: false });
          void navigate({ to: "/new" }).then(() => focusLater(settingFieldId(spec.key)));
        },
      });
    }

    for (const example of catalog?.examples ?? []) {
      const exampleCommand = isAppCommand(example.command) ? example.command : command;

      items.push({
        id: `example-${example.path}`,
        kind: "Example",
        title: example.title,
        description: example.path,
        run: () => {
          loadConfig(example.content, exampleCommand, example.title, false)
            .then((result) => {
              if (!result.loaded) {
                toasts.add({ title: `${example.path} cannot be loaded`, description: result.messages.join("; ") });
              }

              return result;
            })
            .catch((error: Error) => {
              toasts.add({ title: `${example.path} cannot be loaded`, description: error.message });
            });
        },
      });
    }

    return items;
  }, [catalog, command, compareIds, loadConfig, navigate, runList, setTheme, theme, toasts]);
}

function focusLater(id: string) {
  globalThis.setTimeout(() => document.getElementById(id)?.focus(), FOCUS_DELAY_MS);
}
