import { useCallback, useState } from "react";

import { isNumber, isString, type JsonValue } from "../settings/json";
import { parseNumber } from "../settings/numbers";
import { Input } from "../ui/input";

type NumberText = number | string | null;

export function NumberInput({
  value,
  onValueChange,
  ...props
}: Omit<React.ComponentProps<typeof Input>, "type" | "value" | "onChange"> & {
  value: JsonValue;
  onValueChange: (next: NumberText) => void;
}) {
  const current = isNumber(value) || isString(value) ? value : null;
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
  return text.trim() === "" ? null : parseNumber(text);
}

function numberText(value: NumberText): string {
  return value === null ? "" : `${value}`;
}
