import type { ReactNode } from "react";
import Download from "~icons/lucide/download";
import Play from "~icons/lucide/play";
import Settings2 from "~icons/lucide/settings-2";

import { Button } from "../ui/button";
import { Checkbox } from "../ui/checkbox";
import { CopyButton } from "../ui/copy-button";
import { Field, FieldDescription, FieldLabel } from "../ui/field";
import { Input } from "../ui/input";
import { InputGroup, InputGroupAddon, InputGroupInput } from "../ui/input-group";
import { Label } from "../ui/label";
import { NativeSelect, NativeSelectOption } from "../ui/native-select";
import { RadioGroup, RadioGroupItem } from "../ui/radio-group";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "../ui/select";
import { Slider } from "../ui/slider";
import { Switch } from "../ui/switch";
import { Textarea } from "../ui/textarea";
import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";

const MODELS = [
  { value: "gtr", label: "GTR, inferred from the data" },
  { value: "jc69", label: "JC69" },
  { value: "k80", label: "K80" },
  { value: "hky85", label: "HKY85" },
] as const;

const COALESCENT = [
  { value: "none", label: "None" },
  { value: "constant", label: "Constant population size" },
  { value: "skyline", label: "Skyline" },
] as const;

const GRID_POINTS = [400];

const REPRESENTATION = ["sparse"];

const COMMAND = "treetime timetree --tree tree.nwk --alignment aln.fasta.xz --metadata metadata.tsv";

export function ActionDemo() {
  return (
    <div className="flex flex-col gap-3">
      <div className="flex flex-wrap gap-2">
        <Button>
          <Play aria-hidden />
          Run analysis
        </Button>
        <Button variant="outline">Load example</Button>
        <Button variant="ghost">Cancel</Button>
      </div>
      <div className="flex flex-wrap gap-2">
        <Button variant="destructive">Delete run</Button>
        <Button variant="secondary">Secondary</Button>
        <Button variant="link">View log</Button>
        <Button variant="outline" disabled>
          Export
        </Button>
      </div>
      <div className="flex flex-wrap items-center gap-2">
        <Button size="sm">Small</Button>
        <Button size="lg">Large</Button>
        <IconAction label="Settings">
          <Settings2 aria-hidden />
        </IconAction>
        <IconAction label="Download results">
          <Download aria-hidden />
        </IconAction>
        <CopyButton text={COMMAND} label="Copy command" />
      </div>
    </div>
  );
}

export function FieldDemo() {
  return (
    <div className="flex flex-col gap-4">
      <Field>
        <FieldLabel>Run name</FieldLabel>
        <Input defaultValue="flu/h3n2/200 timetree" />
        <FieldDescription>Shown in the run list and in output file names.</FieldDescription>
      </Field>
      <Field>
        <FieldLabel>Clock rate</FieldLabel>
        <InputGroup>
          <InputGroupInput defaultValue="0.003" inputMode="decimal" />
          <InputGroupAddon align="inline-end">subs/site/year</InputGroupAddon>
        </InputGroup>
        <FieldDescription>Leave empty to estimate the rate from the data.</FieldDescription>
      </Field>
      <Field data-invalid>
        <FieldLabel>Maximum iterations</FieldLabel>
        <Input defaultValue="-1" aria-invalid />
        <FieldDescription className="text-destructive">Enter a whole number of 1 or more.</FieldDescription>
      </Field>
    </div>
  );
}

export function SelectDemo() {
  return (
    <div className="flex flex-col gap-4">
      <Field>
        <FieldLabel>Substitution model</FieldLabel>
        <Select items={MODELS} defaultValue="gtr">
          <SelectTrigger className="w-full">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            {MODELS.map((model) => (
              <SelectItem key={model.value} value={model.value}>
                {model.label}
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
      </Field>
      <Field>
        <FieldLabel>Representation, native</FieldLabel>
        <NativeSelect defaultValue="sparse">
          <NativeSelectOption value="sparse">Sparse</NativeSelectOption>
          <NativeSelectOption value="dense">Dense</NativeSelectOption>
        </NativeSelect>
      </Field>
      <Field>
        <FieldLabel>Notes</FieldLabel>
        <Textarea placeholder="What this run is for" />
      </Field>
    </div>
  );
}

export function ChoiceDemo() {
  return (
    <div className="flex flex-col gap-4">
      <div className="flex flex-col gap-2.5">
        <Label>
          <Checkbox defaultChecked />
          Reroot to the best root-to-tip regression
        </Label>
        <Label>
          <Checkbox />
          Resolve polytomies
        </Label>
      </div>
      <RadioGroup defaultValue="constant" aria-label="Coalescent prior">
        {COALESCENT.map((option) => (
          <Label key={option.value}>
            <RadioGroupItem value={option.value} />
            {option.label}
          </Label>
        ))}
      </RadioGroup>
      <Label>
        <Switch defaultChecked />
        Show 90% intervals
      </Label>
      <div className="flex flex-col gap-2">
        <span className="text-sm font-medium">Grid points</span>
        <Slider aria-label="Grid points" defaultValue={GRID_POINTS} min={100} max={1000} step={50} />
      </div>
      <ToggleGroup variant="outline" spacing={0} defaultValue={REPRESENTATION} aria-label="Representation">
        <ToggleGroupItem value="sparse">Sparse</ToggleGroupItem>
        <ToggleGroupItem value="dense">Dense</ToggleGroupItem>
      </ToggleGroup>
    </div>
  );
}

function IconAction({ label, children }: { label: string; children: ReactNode }) {
  return (
    <Tooltip>
      <TooltipTrigger render={<Button variant="ghost" size="icon" aria-label={label} />}>{children}</TooltipTrigger>
      <TooltipContent>{label}</TooltipContent>
    </Tooltip>
  );
}
