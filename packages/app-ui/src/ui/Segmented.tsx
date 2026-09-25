import { useCallback } from "react";

export interface SegmentedOption<Value extends string> {
  value: Value;
  label: string;
}

interface SegmentedProps<Value extends string> {
  options: ReadonlyArray<SegmentedOption<Value>>;
  value: Value;
  onChange: (value: Value) => void;
  label: string;
}

export function Segmented<Value extends string>({ options, value, onChange, label }: SegmentedProps<Value>) {
  return (
    <fieldset className="border-line-strong bg-surface-1 m-0 inline-flex overflow-hidden rounded-md border p-0">
      <legend className="sr-only">{label}</legend>
      {options.map((option) => (
        <SegmentButton key={option.value} option={option} pressed={option.value === value} onChange={onChange} />
      ))}
    </fieldset>
  );
}

function SegmentButton<Value extends string>({
  option,
  pressed,
  onChange,
}: {
  option: SegmentedOption<Value>;
  pressed: boolean;
  onChange: (value: Value) => void;
}) {
  const select = useCallback(() => onChange(option.value), [onChange, option.value]);

  return (
    <button
      type="button"
      aria-pressed={pressed}
      onClick={select}
      className="text-ink-muted aria-pressed:bg-ink aria-pressed:text-surface-1 border-line px-2.5 py-1 text-xs not-first:border-l aria-pressed:font-bold"
    >
      {option.label}
    </button>
  );
}
