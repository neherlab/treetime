import { withListKeys } from "../listKeys";

describe("keys of list items", () => {
  test("repeated texts get distinct keys, and the first keeps its text", () => {
    const keys = withListKeys(["a.nwk", "b.nwk", "a.nwk", "a.nwk"], (name) => name).map(({ key }) => key);

    expect([new Set(keys).size, keys[0], keys[1]]).toStrictEqual([4, "a.nwk", "b.nwk"]);
  });
});
