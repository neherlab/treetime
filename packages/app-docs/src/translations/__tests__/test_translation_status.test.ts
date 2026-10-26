import { contentHash, translationStatus } from "../translation-status";

const ENGLISH_PAGE = "---\ntitle: Setup\n---\n\nInstall the tools.\n";

const ENGLISH_META = '{ "title": "Development", "pages": ["setup"] }\n';

describe("translation status", () => {
  test("a translation that records the hash of its English page is current", () => {
    const files = new Map([
      ["en/setup.md", ENGLISH_PAGE],
      ["de/setup.md", `---\ntitle: Einrichtung\nsourceHash: ${contentHash(ENGLISH_PAGE)}\n---\n\nText.\n`],
    ]);

    expect(translationStatus(files)).toStrictEqual([]);
  });

  test("lists a missing translation with the hash to record", () => {
    const files = new Map([["en/setup.md", ENGLISH_PAGE]]);

    expect(translationStatus(files)).toStrictEqual([
      { status: "missing", path: "de/setup.md", sourceHash: contentHash(ENGLISH_PAGE) },
    ]);
  });

  test("lists a translation of an older English page, or without a hash, as stale", () => {
    const files = new Map([
      ["en/setup.md", ENGLISH_PAGE],
      ["en/meta.json", ENGLISH_META],
      ["de/setup.md", "---\ntitle: Einrichtung\nsourceHash: 0123456789abcdef\n---\n\nText.\n"],
      ["de/meta.json", '{ "title": "Entwicklung", "pages": ["setup"] }\n'],
    ]);

    expect(translationStatus(files).map(({ status, path }) => `${status} ${path}`)).toStrictEqual([
      "stale de/meta.json",
      "stale de/setup.md",
    ]);
  });

  test("reads the hash of a navigation file from its JSON", () => {
    const files = new Map([
      ["en/meta.json", ENGLISH_META],
      ["de/meta.json", JSON.stringify({ title: "Entwicklung", sourceHash: contentHash(ENGLISH_META) })],
    ]);

    expect(translationStatus(files)).toStrictEqual([]);
  });

  test("lists a translation without an English page as an orphan", () => {
    const files = new Map([["de/removed.md", "---\ntitle: Alt\n---\n"]]);

    expect(translationStatus(files)).toStrictEqual([
      { status: "orphan", path: "de/removed.md", sourceHash: undefined },
    ]);
  });
});
