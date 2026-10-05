import type { JsonValue } from "@neherlab/app-contracts";
import { useCallback, useState } from "react";

import { isNumber, isString } from "../settings/json";
import { parseNumber } from "../settings/numbers";
import { Input } from "../ui/input";

type NumberText = number | string | undefined;

export function NumberInput({
  value,
  onValueChange,
  ...props
}: Omit<React.ComponentProps<typeof Input>, "type" | "value" | "onChange"> & {
  value: JsonValue | undefined;
  onValueChange: (next: NumberText) => void;
}) {
  const current = isNumber(value) || isString(value) ? value : undefined;
  const [text, setText] = useState(() => numberText(current));
  const shown = parse(text) === current ? text : numberText(current);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      setText(event.target.value);
      onValueChange(parse(event.target.value));
    },
    [onValueChange],
  );

  return <Input type="text" autoComplete="off" spellCheck={false} value={shown} onChange={onChange} {...props} />;
}

function parse(text: string): NumberText {
  return text.trim() === "" ? undefined : parseNumber(text);
}

function numberText(value: NumberText): string {
  return value === undefined ? "" : `${value}`;
}
