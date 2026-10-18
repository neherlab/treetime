import { filesForSetting, storedSource, storedValue, storeInOrder } from "../fileSelection";

describe("files chosen for an input setting", () => {
  test("a list setting keeps every file and a single setting keeps the first", () => {
    expect([filesForSetting(["a", "b"], true), filesForSetting(["a", "b"], false)]).toStrictEqual([["a", "b"], ["a"]]);
  });

  test("files are stored one at a time in selection order", async () => {
    const events: string[] = [];

    const stored = await storeInOrder(
      ["ha", "na"],
      async (name) => {
        events.push(`start ${name}`);
        await Promise.resolve();
        events.push(`end ${name}`);

        return { name, path: `/inputs/${name}`, size: name.length };
      },
      () => {
        events.push("discard");
      },
    );

    expect([events, stored.map((file) => file.path)]).toStrictEqual([
      ["start ha", "end ha", "start na", "end na"],
      ["/inputs/ha", "/inputs/na"],
    ]);
  });

  test("the setting value and the source label follow the stored files", () => {
    const stored = [
      { name: "aligned.fasta", path: "/inputs/1/aligned.fasta", size: 8 },
      { name: "aligned.fasta", path: "/inputs/2/aligned.fasta", size: 6 },
    ];

    expect([
      storedValue(stored, true),
      storedValue(stored, false),
      storedValue([], true),
      storedSource(stored),
    ]).toStrictEqual([
      ["/inputs/1/aligned.fasta", "/inputs/2/aligned.fasta"],
      "/inputs/1/aligned.fasta",
      undefined,
      { label: "aligned.fasta, aligned.fasta", size: 14 },
    ]);
  });
});
