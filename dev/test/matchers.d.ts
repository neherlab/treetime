declare module "bun:test" {
  interface Matchers<T> {
    // oxlint-disable-next-line anti-slop/no-unknown-parameters -- an equality matcher compares against a value of any type, as in vitest
    toEqual(expected: unknown): void;
    // oxlint-disable-next-line anti-slop/no-unknown-parameters -- an equality matcher compares against a value of any type, as in vitest
    toStrictEqual(expected: unknown): void;
  }
}
