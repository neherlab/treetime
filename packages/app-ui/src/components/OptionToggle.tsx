import { useCallback, useMemo } from "react";

import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { toggledChoice } from "./toggleChoice";

export interface ToggleOption<Value extends string> {
  value: Value;
  label: string;
}

export function OptionToggle<Value extends string>({
  options,
  value,
  onChange,
  label,
}: {
  options: ReadonlyArray<ToggleOption<Value>>;
  value: Value;
  onChange: (value: Value) => void;
  label: string;
}) {
  const pressed = useMemo(() => [value], [value]);
  const values = useMemo(() => options.map((option) => option.value), [options]);

  const onValueChange = useCallback(
    (next: string[]) => {
      const picked = toggledChoice(values, value, next);

      if (picked !== undefined) {
        onChange(picked);
      }
    },
    [onChange, value, values],
  );

  return (
    <ToggleGroup
      aria-label={label}
      variant="outline"
      size="sm"
      spacing={0}
      value={pressed}
      onValueChange={onValueChange}
    >
      {options.map((option) => (
        <ToggleGroupItem key={option.value} value={option.value}>
          {option.label}
        </ToggleGroupItem>
      ))}
    </ToggleGroup>
  );
}
