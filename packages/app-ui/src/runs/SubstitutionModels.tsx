import type {
  AlignmentSummary,
  RunModel,
  RunModels,
  SequenceModelSummary,
  TraitModelSummary,
} from "@neherlab/app-contracts";
import { useMemo } from "react";

import { Panel } from "../components/Panel";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";

export function SubstitutionModels({ models: { alignments, models } }: { models: RunModels }) {
  const sequences = useMemo(
    () => models.filter((model): model is RunModel & { kind: "sequence" } => model.kind === "sequence"),
    [models],
  );

  const traits = useMemo(
    () => models.filter((model): model is RunModel & { kind: "trait" } => model.kind === "trait"),
    [models],
  );

  return (
    <Panel title="Substitution models" hint="The models of the run's model file">
      {sequences.length > 0 && <SequenceModels models={sequences} />}
      {sequences.length > 0 && alignments.length > 0 && <Alignments alignments={alignments} />}
      {traits.length > 0 && <TraitModels models={traits} />}
    </Panel>
  );
}

function SequenceModels({ models }: { models: readonly SequenceModelSummary[] }) {
  return (
    <Table aria-label="Sequence models">
      <TableHeader>
        <TableRow>
          <TableHead className="text-right">Partition</TableHead>
          <TableHead>Name</TableHead>
          <TableHead className="text-right">Alignment</TableHead>
          <TableHead>File</TableHead>
          <TableHead className="text-right">Sites</TableHead>
          <TableHead>Model</TableHead>
          <TableHead className="text-right">Overall rate μ</TableHead>
          <TableHead>States, in the order of W and π</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {models.map((model) => (
          <TableRow key={model.partition}>
            <TableCell className="text-right">{model.partition}</TableCell>
            <TableCell>{model.name ?? "-"}</TableCell>
            <TableCell className="text-right">{model.alignment}</TableCell>
            <TableCell>
              <code className="font-mono text-xs">{model.path ?? "-"}</code>
            </TableCell>
            <TableCell className="text-right">{model.sites}</TableCell>
            <TableCell>{model.model}</TableCell>
            <TableCell className="text-right">{model.mu.toPrecision(4)}</TableCell>
            <TableCell className="text-xs whitespace-normal">{model.states.join(", ")}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}

function Alignments({ alignments }: { alignments: readonly AlignmentSummary[] }) {
  return (
    <Table aria-label="Alignment files">
      <TableHeader>
        <TableRow>
          <TableHead className="text-right">Alignment</TableHead>
          <TableHead>File</TableHead>
          <TableHead className="text-right">Sites</TableHead>
          <TableHead className="text-right">Indel rate</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {alignments.map((alignment) => (
          <TableRow key={alignment.alignment}>
            <TableCell className="text-right">{alignment.alignment}</TableCell>
            <TableCell>
              <code className="font-mono text-xs">{alignment.path}</code>
            </TableCell>
            <TableCell className="text-right">{alignment.sites}</TableCell>
            <TableCell className="text-right">{alignment.indel_rate?.toPrecision(4) ?? "-"}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}

function TraitModels({ models }: { models: readonly TraitModelSummary[] }) {
  return (
    <Table aria-label="Trait models">
      <TableHeader>
        <TableRow>
          <TableHead>Attribute</TableHead>
          <TableHead>Model</TableHead>
          <TableHead className="text-right">Overall rate μ</TableHead>
          <TableHead>States, in the order of W and π</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {models.map((model) => (
          <TableRow key={model.attribute}>
            <TableCell>{model.attribute}</TableCell>
            <TableCell>{model.model}</TableCell>
            <TableCell className="text-right">{model.mu.toPrecision(4)}</TableCell>
            <TableCell className="text-xs whitespace-normal">{model.states.join(", ")}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}
