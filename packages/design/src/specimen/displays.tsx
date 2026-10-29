import CircleCheck from "~icons/lucide/circle-check";
import CircleX from "~icons/lucide/circle-x";
import Ellipsis from "~icons/lucide/ellipsis";
import Info from "~icons/lucide/info";
import TriangleAlert from "~icons/lucide/triangle-alert";

import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Badge } from "../ui/badge";
import {
  Breadcrumb,
  BreadcrumbItem,
  BreadcrumbLink,
  BreadcrumbList,
  BreadcrumbPage,
  BreadcrumbSeparator,
} from "../ui/breadcrumb";
import { Button } from "../ui/button";
import { DescriptionDetails, DescriptionList, DescriptionTerm } from "../ui/description-list";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle, DialogTrigger } from "../ui/dialog";
import { DropdownMenu, DropdownMenuContent, DropdownMenuItem, DropdownMenuTrigger } from "../ui/dropdown-menu";
import { Kbd } from "../ui/kbd";
import { Meter, MeterLabel, MeterValue } from "../ui/meter";
import {
  Pagination,
  PaginationContent,
  PaginationItem,
  PaginationLink,
  PaginationNext,
  PaginationPrevious,
} from "../ui/pagination";
import {
  Popover,
  PopoverContent,
  PopoverDescription,
  PopoverHeader,
  PopoverTitle,
  PopoverTrigger,
} from "../ui/popover";
import { Progress, ProgressLabel, ProgressValue } from "../ui/progress";
import { Skeleton } from "../ui/skeleton";
import { Spinner } from "../ui/spinner";
import { Stat } from "../ui/stat";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { Tabs, TabsList, TabsTrigger } from "../ui/tabs";

const NODES = [
  { name: "A/Texas/50/2019", date: "2019.88", interval: "sampled", mutations: 2, clade: "3C.2a1b.1b" },
  { name: "NODE_0000142", date: "2019.31", interval: "2019.05 to 2019.58", mutations: 4, clade: "3C.2a1b.1b" },
  { name: "A/Hong_Kong/4801/2018", date: "2018.97", interval: "sampled", mutations: 1, clade: "3C.3a" },
  { name: "NODE_0000087", date: "2016.42", interval: "2016.18 to 2016.61", mutations: 3, clade: "3C.2a" },
] as const;

const SELECTED_NODE = "NODE_0000142";

export function StateDemo() {
  return (
    <div className="flex flex-col gap-4">
      <Tabs defaultValue="results">
        <TabsList variant="line">
          <TabsTrigger value="results">Results</TabsTrigger>
          <TabsTrigger value="log">Log</TabsTrigger>
          <TabsTrigger value="settings">Settings</TabsTrigger>
        </TabsList>
      </Tabs>
      <div className="flex flex-wrap gap-2">
        <Badge variant="outline">v1.0.0-nightly</Badge>
        <Badge>Sparse</Badge>
        <Badge variant="success">
          <CircleCheck aria-hidden />
          Converged
        </Badge>
        <Badge variant="warning">
          <TriangleAlert aria-hidden />3 warnings
        </Badge>
        <Badge variant="destructive">
          <CircleX aria-hidden />
          Failed
        </Badge>
      </div>
      <ul className="flex flex-col gap-1.5 text-sm">
        <li className="flex items-center gap-2">
          <Spinner className="text-primary size-3.5" aria-hidden />
          Running, iteration 2 of 3
        </li>
        <li className="flex items-center gap-2">
          <CircleCheck className="text-success size-3.5" aria-hidden />
          Finished in 14 s
        </li>
        <li className="flex items-center gap-2">
          <CircleX className="text-destructive size-3.5" aria-hidden />
          Failed: the tree file cannot be read
        </li>
      </ul>
    </div>
  );
}

export function NoticeDemo() {
  return (
    <div className="flex flex-col gap-3">
      <Alert>
        <Info aria-hidden />
        <AlertTitle>Rate from root-to-tip regression</AlertTitle>
        <AlertDescription>Set a clock rate to skip the estimate.</AlertDescription>
      </Alert>
      <Alert variant="success">
        <CircleCheck aria-hidden />
        <AlertTitle>Run finished</AlertTitle>
        <AlertDescription>12 output files in out/flu-h3n2-200.</AlertDescription>
      </Alert>
      <Alert variant="warning">
        <TriangleAlert aria-hidden />
        <AlertTitle>3 sequences have no date</AlertTitle>
        <AlertDescription>They stay in the tree but not in the clock fit.</AlertDescription>
      </Alert>
      <Alert variant="destructive">
        <CircleX aria-hidden />
        <AlertTitle>The tree file cannot be read</AlertTitle>
        <AlertDescription>Line 1, column 214: unexpected ). Check that the parentheses balance.</AlertDescription>
      </Alert>
    </div>
  );
}

export function ProgressDemo() {
  return (
    <div className="flex flex-col gap-4">
      <Progress value={66}>
        <ProgressLabel>Optimizing branch lengths</ProgressLabel>
        <ProgressValue />
      </Progress>
      <Meter value={15}>
        <MeterLabel>Memory</MeterLabel>
        <MeterValue />
      </Meter>
      <div className="flex items-center gap-3">
        <Skeleton className="size-9 rounded-full" />
        <div className="flex flex-1 flex-col gap-2">
          <Skeleton className="h-3 w-3/4" />
          <Skeleton className="h-3 w-1/2" />
        </div>
      </div>
    </div>
  );
}

export function ReadoutDemo() {
  return (
    <div className="flex flex-col gap-4">
      <div className="grid grid-cols-2 gap-4">
        <Stat label="Clock rate" value="3.54" uncertainty="± 0.21 × 10⁻³" detail="subs/site/year" />
        <Stat label="Root date" value="2009.4" detail="90%: 2008.9 to 2009.8" />
      </div>
      <DescriptionList>
        <DescriptionTerm>Tree</DescriptionTerm>
        <DescriptionDetails>
          <code className="font-mono">tree.nwk</code>, 200 tips
        </DescriptionDetails>
        <DescriptionTerm>Alignment</DescriptionTerm>
        <DescriptionDetails>
          <code className="font-mono">aln.fasta.xz</code>, 1,701 sites
        </DescriptionDetails>
        <DescriptionTerm>Model</DescriptionTerm>
        <DescriptionDetails>GTR, inferred</DescriptionDetails>
      </DescriptionList>
      <p className="text-muted-foreground flex items-center gap-1.5 text-sm">
        Command palette <Kbd>Ctrl</Kbd>
        <Kbd>K</Kbd>
      </p>
    </div>
  );
}

export function OverlayDemo() {
  return (
    <div className="flex flex-wrap gap-2">
      <DropdownMenu>
        <DropdownMenuTrigger render={<Button variant="outline" aria-label="Run actions" />}>
          <Ellipsis aria-hidden />
          Actions
        </DropdownMenuTrigger>
        <DropdownMenuContent>
          <DropdownMenuItem>Open output folder</DropdownMenuItem>
          <DropdownMenuItem>Copy command</DropdownMenuItem>
          <DropdownMenuItem>Run again</DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
      <Popover>
        <PopoverTrigger render={<Button variant="outline" />}>Clock rate</PopoverTrigger>
        <PopoverContent>
          <PopoverHeader>
            <PopoverTitle>Clock rate</PopoverTitle>
            <PopoverDescription>Substitutions per site per year, from the root-to-tip regression.</PopoverDescription>
          </PopoverHeader>
        </PopoverContent>
      </Popover>
      <Dialog>
        <DialogTrigger render={<Button variant="destructive" />}>Delete run...</DialogTrigger>
        <DialogContent>
          <DialogHeader>
            <DialogTitle>Delete this run?</DialogTitle>
            <DialogDescription>
              The run and its 12 output files are removed from this computer. You cannot undo this.
            </DialogDescription>
          </DialogHeader>
        </DialogContent>
      </Dialog>
    </div>
  );
}

export function NavigationDemo() {
  return (
    <div className="flex flex-col gap-4">
      <Breadcrumb>
        <BreadcrumbList>
          <BreadcrumbItem>
            <BreadcrumbLink href="#runs">Runs</BreadcrumbLink>
          </BreadcrumbItem>
          <BreadcrumbSeparator />
          <BreadcrumbItem>
            <BreadcrumbLink href="#run">flu/h3n2/200 timetree</BreadcrumbLink>
          </BreadcrumbItem>
          <BreadcrumbSeparator />
          <BreadcrumbItem>
            <BreadcrumbPage>Nodes</BreadcrumbPage>
          </BreadcrumbItem>
        </BreadcrumbList>
      </Breadcrumb>
      <Pagination className="justify-start">
        <PaginationContent>
          <PaginationItem>
            <PaginationPrevious href="#previous" />
          </PaginationItem>
          <PaginationItem>
            <PaginationLink href="#page-1" isActive>
              1
            </PaginationLink>
          </PaginationItem>
          <PaginationItem>
            <PaginationLink href="#page-2">2</PaginationLink>
          </PaginationItem>
          <PaginationItem>
            <PaginationNext href="#next" />
          </PaginationItem>
        </PaginationContent>
      </Pagination>
    </div>
  );
}

export function TableDemo() {
  return (
    <div className="overflow-x-auto rounded-lg border">
      <Table>
        <TableHeader>
          <TableRow>
            <TableHead>Node</TableHead>
            <TableHead className="text-right">Date</TableHead>
            <TableHead className="text-right">90% interval</TableHead>
            <TableHead className="text-right">Mutations</TableHead>
            <TableHead>Clade</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {NODES.map((node) => (
            <TableRow key={node.name} data-state={node.name === SELECTED_NODE ? "selected" : undefined}>
              <TableCell className="font-mono text-xs">{node.name}</TableCell>
              <TableCell className="text-right tabular-nums">{node.date}</TableCell>
              <TableCell className="text-muted-foreground text-right tabular-nums">{node.interval}</TableCell>
              <TableCell className="text-right tabular-nums">{node.mutations}</TableCell>
              <TableCell>{node.clade}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
}
