import type { RunModel, SequenceModelSummary, TraitModelSummary } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { Panel } from "../components/Panel";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";

export function SubstitutionModels({ models }: { models: readonly RunModel[] }) {
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
      {traits.length > 0 && <TraitModels models={traits} />}
    </Panel>
  );
}

function SequenceModels({ models }: { models: readonly SequenceModelSummary[] }) {
  return (
    <Table aria-label="Sequence models">
      <TableHeader>
        <TableRow>
          <TableHead className="text-right">Alignment</TableHead>
          <TableHead>File</TableHead>
          <TableHead>Model</TableHead>
          <TableHead className="text-right">Overall rate μ</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {models.map((model) => (
          <TableRow key={model.partition}>
            <TableCell className="text-right">{model.partition}</TableCell>
            <TableCell>
              <code className="font-mono text-xs">{model.alignment}</code>
            </TableCell>
            <TableCell>{model.name}</TableCell>
            <TableCell className="text-right">{model.mu.toPrecision(4)}</TableCell>
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
            <TableCell>{model.name}</TableCell>
            <TableCell className="text-right">{model.mu.toPrecision(4)}</TableCell>
            <TableCell className="text-xs whitespace-normal">{model.states.join(", ")}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}
