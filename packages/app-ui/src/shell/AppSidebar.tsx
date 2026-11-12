import type { AppCommand, RunSummary } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Checkbox } from "@neherlab/design/ui/checkbox";
import { DropdownMenuItem } from "@neherlab/design/ui/dropdown-menu";
import { InputGroup, InputGroupAddon, InputGroupInput } from "@neherlab/design/ui/input-group";
import { Kbd } from "@neherlab/design/ui/kbd";
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarGroup,
  SidebarHeader,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
} from "@neherlab/design/ui/sidebar";
import { ToggleGroup, ToggleGroupItem } from "@neherlab/design/ui/toggle-group";
import { useHotkey } from "@tanstack/react-hotkeys";
import { Link, useParams } from "@tanstack/react-router";
import { cn } from "cn";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";
import Plus from "~icons/lucide/plus";
import Search from "~icons/lucide/search";

import { useNewAnalysis } from "../analysis/useDraftSession";
import { toggledChoice } from "../components/toggleChoice";
import { dayLabel, formatBytes, headlineText } from "../format";
import { COMMANDS, commandSettings } from "../settings/catalog";
import { useShellStore } from "../store/shell";
import { RunActionsMenu, useRunActions } from "./runActions";
import { groupRuns, type RunEntry } from "./runList";
import { StatusIcon, statusLabel } from "./StatusIcon";
import { RUN_FILTER_ID } from "./useGlobalShortcuts";
import { useListedRuns } from "./useListedRuns";

const ALL_COMMANDS = "all";

const DEPTH_INDENT = ["", "pl-6", "pl-10", "pl-14"] as const;

export function AppSidebar() {
  const runFilter = useShellStore((state) => state.runFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const { all, runs, loaded, error } = useListedRuns(runFilter, commandFilter);
  const groups = useMemo(() => groupRuns(runs), [runs]);

  const commands = useMemo(
    () => COMMANDS.flatMap(({ command }) => (all.some((run) => run.command === command) ? [command] : [])),
    [all],
  );

  return (
    <Sidebar className="top-(--header-height) h-[calc(100svh-var(--header-height))]!">
      <SidebarHeader className="gap-2">
        <NewAnalysisButton />
        <RunFilter />
        {commands.length > 1 && <CommandFilter commands={commands} />}
      </SidebarHeader>
      <SidebarContent className="overscroll-contain">
        {error !== null && (
          <p role="alert" className="text-destructive px-4 py-2 text-xs">
            The runs cannot be listed: {error.message}
          </p>
        )}
        <nav aria-label="Runs" className="flex flex-col">
          {groups.map((group) => (
            <SidebarGroup key={group.id} className="border-b py-1.5 last:border-b-0">
              <SidebarMenu aria-label={group.entries[0]?.run.title}>
                {group.entries.map((entry) => (
                  <RunItem key={entry.run.id} entry={entry} />
                ))}
              </SidebarMenu>
            </SidebarGroup>
          ))}
        </nav>
        {groups.length === 0 && loaded && (
          <p className="text-muted-foreground px-4 py-6 text-center text-sm">
            {all.length === 0 ? "No runs yet. Start a new analysis." : "No run matches the filter."}
          </p>
        )}
      </SidebarContent>
      <SelectionFooter runs={all} />
    </Sidebar>
  );
}

function NewAnalysisButton() {
  const openNew = useNewAnalysis();

  return (
    <Button className="w-full" onClick={openNew}>
      <Plus aria-hidden />
      New analysis
      <Kbd className="bg-primary-foreground/15 text-primary-foreground ml-auto">N</Kbd>
    </Button>
  );
}

function RunFilter() {
  const runFilter = useShellStore((state) => state.runFilter);
  const setRunFilter = useShellStore((state) => state.setRunFilter);

  const onFilter = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => setRunFilter(event.target.value),
    [setRunFilter],
  );

  return (
    <InputGroup className="bg-background h-8">
      <InputGroupAddon>
        <Search aria-hidden />
      </InputGroupAddon>
      <InputGroupInput
        id={RUN_FILTER_ID}
        type="search"
        value={runFilter}
        onChange={onFilter}
        placeholder="Filter by name"
        aria-label="Filter runs"
      />
      <InputGroupAddon align="inline-end">
        <Kbd>/</Kbd>
      </InputGroupAddon>
    </InputGroup>
  );
}

function CommandFilter({ commands }: { commands: readonly AppCommand[] }) {
  const commandFilter = useShellStore((state) => state.commandFilter);
  const setCommandFilter = useShellStore((state) => state.setCommandFilter);
  const pressed = useMemo(() => [commandFilter ?? ALL_COMMANDS], [commandFilter]);

  const onValueChange = useCallback(
    (next: string[]) => {
      const picked = toggledChoice<string>([ALL_COMMANDS, ...commands], commandFilter ?? ALL_COMMANDS, next);

      if (picked !== undefined) {
        setCommandFilter(commands.find((command) => command === picked) ?? null);
      }
    },
    [commandFilter, commands, setCommandFilter],
  );

  return (
    <ToggleGroup
      aria-label="Filter by analysis"
      size="sm"
      spacing={1}
      value={pressed}
      onValueChange={onValueChange}
      className="flex-wrap"
    >
      <ToggleGroupItem value={ALL_COMMANDS} className="h-6 px-2 text-xs">
        All
      </ToggleGroupItem>
      {commands.map((command) => (
        <ToggleGroupItem key={command} value={command} className="h-6 px-2 text-xs">
          {commandSettings(command).title}
        </ToggleGroupItem>
      ))}
    </ToggleGroup>
  );
}

function RunItem({ entry: { run, depth } }: { entry: RunEntry }) {
  const currentId = useParams({ strict: false, select: (params) => params.id });
  const selecting = useShellStore((state) => state.selectedIds.length > 0);
  const selected = useShellStore((state) => state.selectedIds.includes(run.id));
  const toggleSelected = useShellStore((state) => state.toggleSelected);
  const headline = run.status === "ok" ? headlineText(run.headline) : "";
  const onToggle = useCallback(() => toggleSelected(run.id), [run.id, toggleSelected]);

  const onClick = useCallback(
    (event: React.MouseEvent) => {
      if (selecting || event.metaKey || event.ctrlKey || event.shiftKey) {
        event.preventDefault();
        toggleSelected(run.id);
      }
    },
    [run.id, selecting, toggleSelected],
  );

  return (
    <SidebarMenuItem className={DEPTH_INDENT[Math.min(depth, DEPTH_INDENT.length - 1)]}>
      <SidebarMenuButton
        isActive={run.id === currentId}
        className={cn("h-auto items-start py-1.5 pr-8", selected && "bg-sidebar-accent ring-primary/40 ring-1")}
        render={<Link to="/runs/$id/results" params={{ id: run.id }} onClick={onClick} />}
      >
        {selecting ? (
          <Checkbox
            checked={selected}
            onCheckedChange={onToggle}
            aria-label={`Select ${run.title}`}
            className="mt-0.5"
          />
        ) : (
          <span className="mt-0.5" title={statusLabel(run.status)}>
            <StatusIcon status={run.status} />
          </span>
        )}
        <span className="grid min-w-0 flex-1 gap-0.5">
          <span className="flex items-baseline gap-2">
            <span className="truncate font-bold" title={run.title}>
              {run.title}
            </span>
            <span className="text-muted-foreground ml-auto text-xs whitespace-nowrap">{headline}</span>
          </span>
          <span className="text-muted-foreground truncate text-xs">
            {dayLabel(run.created_at, DateTime.now())}, {formatBytes(run.size)}
          </span>
        </span>
      </SidebarMenuButton>
      <div className="absolute top-1.5 right-1">
        <RunActionsMenu run={run}>
          <DropdownMenuItem onClick={onToggle}>{selected ? "Deselect" : "Select"}</DropdownMenuItem>
        </RunActionsMenu>
      </div>
    </SidebarMenuItem>
  );
}

function SelectionFooter({ runs }: { runs: readonly RunSummary[] }) {
  const selectedIds = useShellStore((state) => state.selectedIds);
  const clearSelection = useShellStore((state) => state.clearSelection);
  const actions = useRunActions();
  const selected = useMemo(() => runs.filter((run) => selectedIds.includes(run.id)), [runs, selectedIds]);
  const [first, second] = selected;
  const comparable = selected.length === 2 && selected.every((run) => run.status === "ok");
  const onDelete = useCallback(() => actions.remove(selected), [actions, selected]);

  useHotkey("Escape", clearSelection, { enabled: selected.length > 0 });

  if (selected.length === 0) {
    return null;
  }

  return (
    <SidebarFooter className="border-t">
      <div className="flex items-center gap-2 text-xs">
        <span className="text-muted-foreground flex-1" aria-live="polite">
          {selected.length === 1 ? "1 run selected" : `${selected.length} runs selected`}
        </span>
        {comparable && first !== undefined && second !== undefined && (
          <Button
            size="xs"
            nativeButton={false}
            render={<Link to="/compare/$a/$b" params={{ a: first.id, b: second.id }} />}
          >
            Compare
          </Button>
        )}
        <Button size="xs" variant="outline" onClick={onDelete}>
          Delete
        </Button>
        <Button size="xs" variant="ghost" onClick={clearSelection}>
          Clear
        </Button>
      </div>
      {selected.length === 1 && (
        <p className="text-muted-foreground text-xs">Select one more finished run to compare.</p>
      )}
    </SidebarFooter>
  );
}
