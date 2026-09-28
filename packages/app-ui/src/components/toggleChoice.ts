export function toggledChoice<Value extends string>(
  values: readonly Value[],
  current: Value,
  pressed: readonly string[],
): Value | undefined {
  return values.find((value) => value !== current && pressed.includes(value));
}
