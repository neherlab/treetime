import type { AppCommand, RunSummary } from "@neherlab/app-contracts";
import { runsList } from "@neherlab/app-contracts/client";
import { Link, useNavigate, useParams } from "@tanstack/react-router";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";
import Plus from "~icons/lucide/plus";
import Search from "~icons/lucide/search";

import { useApi } from "../api/hooks";
import { toggledChoice } from "../components/toggleChoice";
import { headlineText } from "../format";
import { COMMANDS, commandSettings } from "../settings/catalog";
import { useShellStore } from "../store/shell";
import { Badge } from "../ui/badge";
import { Button } from "../ui/button";
import { Checkbox } from "../ui/checkbox";
import { InputGroup, InputGroupAddon, InputGroupInput } from "../ui/input-group";
import { Kbd } from "../ui/kbd";
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarGroup,
  SidebarGroupLabel,
  SidebarHeader,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
} from "../ui/sidebar";
import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { changedFlags, groupRuns, listedRuns } from "./runList";
import { StatusIcon, statusLabel } from "./StatusIcon";
import { RUN_FILTER_ID } from "./useGlobalShortcuts";

const FLAGS_SHOWN = 2;

const ALL_COMMANDS = "all";

export function AppSidebar() {
  const { data, error } = useApi((context) => runsList(context));
  const runFilter = useShellStore((state) => state.runFilter);
  const commandFilter = useShellStore((state) => state.commandFilter);
  const allRuns = useMemo(() => listedRuns(data?.runs ?? [], "", null), [data]);
  const runs = useMemo(() => listedRuns(allRuns, runFilter, commandFilter), [allRuns, runFilter, commandFilter]);
  const groups = useMemo(() => groupRuns(runs, DateTime.now()), [runs]);

  const commands = useMemo(
    () => COMMANDS.flatMap(({ command }) => (allRuns.some((run) => run.command === command) ? [command] : [])),
    [allRuns],
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
        <nav aria-label="Runs" className="flex flex-col gap-2">
          {groups.map((group) => (
            <SidebarGroup key={group.label}>
              <SidebarGroupLabel render={renderHeading} className="justify-between">
                <span>{group.label}</span>
                <span className="tabular-nums">{group.runs.length}</span>
              </SidebarGroupLabel>
              <SidebarMenu>
                {group.runs.map((run) => (
                  <RunItem key={run.id} run={run} />
                ))}
              </SidebarMenu>
            </SidebarGroup>
          ))}
        </nav>
        {groups.length === 0 && data !== undefined && (
          <p className="text-muted-foreground px-4 py-6 text-center text-sm">
            {allRuns.length === 0 ? "No runs yet. Start a new analysis." : "No run matches the filter."}
          </p>
        )}
      </SidebarContent>
      <CompareFooter />
    </Sidebar>
  );
}

function NewAnalysisButton() {
  const navigate = useNavigate();
  const openNew = useCallback(() => void navigate({ to: "/new" }), [navigate]);

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
        placeholder="Filter by name or setting"
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

function RunItem({ run }: { run: RunSummary }) {
  const currentId = useParams({ strict: false, select: (params) => params.id });
  const flags = changedFlags(run);
  const headline = run.status === "ok" ? headlineText(run.headline) : "";

  return (
    <SidebarMenuItem>
      <SidebarMenuButton
        isActive={run.id === currentId}
        className="h-auto items-start py-1.5 pr-8"
        render={<Link to="/runs/$id/results" params={{ id: run.id }} />}
      >
        <span className="mt-0.5" title={statusLabel(run.status)}>
          <StatusIcon status={run.status} />
        </span>
        <span className="grid min-w-0 flex-1 gap-0.5">
          <span className="flex items-baseline gap-2">
            <span className="truncate font-bold" title={run.title}>
              {run.title}
            </span>
            <span className="text-muted-foreground ml-auto text-xs whitespace-nowrap">{headline}</span>
          </span>
          <span className="text-muted-foreground flex items-center gap-1 overflow-hidden text-xs">
            <span className="text-primary font-bold">{run.command}</span>
            {flags.slice(0, FLAGS_SHOWN).map((flag) => (
              <Badge
                key={flag}
                variant="secondary"
                title={flag}
                className="text-muted-foreground h-4 max-w-[26ch] rounded-sm px-1 font-mono font-normal"
              >
                <span className="min-w-0 truncate">{flag}</span>
              </Badge>
            ))}
            {flags.length > FLAGS_SHOWN && <span>+{flags.length - FLAGS_SHOWN}</span>}
          </span>
        </span>
      </SidebarMenuButton>
      {run.status === "ok" && <CompareToggle run={run} />}
    </SidebarMenuItem>
  );
}

function CompareToggle({ run }: { run: RunSummary }) {
  const compareIds = useShellStore((state) => state.compareIds);
  const toggleCompare = useShellStore((state) => state.toggleCompare);
  const onCompare = useCallback(() => toggleCompare(run.id), [run.id, toggleCompare]);

  return (
    <Checkbox
      checked={compareIds.includes(run.id)}
      onCheckedChange={onCompare}
      title="Select to compare"
      aria-label={`Select ${run.title} to compare`}
      className="bg-sidebar absolute right-2 bottom-2"
    />
  );
}

function CompareFooter() {
  const compareIds = useShellStore((state) => state.compareIds);
  const clearCompare = useShellStore((state) => state.clearCompare);
  const [first, second] = compareIds;

  if (first === undefined) {
    return null;
  }

  return (
    <SidebarFooter className="border-t">
      <div className="flex items-center gap-2 text-xs">
        <span className="text-muted-foreground flex-1">
          {second === undefined ? "Select one more run to compare" : "2 runs selected"}
        </span>
        {second !== undefined && (
          <Button size="xs" nativeButton={false} render={<Link to="/compare/$a/$b" params={{ a: first, b: second }} />}>
            Compare
          </Button>
        )}
        <Button size="xs" variant="ghost" onClick={clearCompare}>
          Clear
        </Button>
      </div>
    </SidebarFooter>
  );
}

function renderHeading({ children, ...props }: React.ComponentProps<"h2">) {
  return <h2 {...props}>{children}</h2>;
}
