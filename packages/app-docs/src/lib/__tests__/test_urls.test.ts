import { assetPath, docsPath } from "../urls";
import { versionPageUrl } from "../version-links";

describe("docs urls", () => {
  test.each([
    ["the English home page", [], "en", "/dev"],
    ["an English page", ["showcase", "math"], "en", "/dev/showcase/math"],
    ["the German home page", [], "de", "/dev/de"],
    ["a German page", ["showcase", "math"], "de", "/dev/de/showcase/math"],
  ] as const)("%s", (_, slugs, locale, expected) => {
    expect(docsPath("/dev", slugs, locale)).toBe(expected);
  });

  test("a file of the public folder lies below the version", () => {
    expect([assetPath("/dev", "favicon.svg"), assetPath("/dev", "/favicon.svg")]).toStrictEqual([
      "/dev/favicon.svg",
      "/dev/favicon.svg",
    ]);
  });
});

describe("version page urls", () => {
  const manifest = { version: "1.0", pages: { en: ["", "showcase/math"], de: [""] } };

  test.each([
    ["a page the target version has", "en", "showcase/math", "/1.0/showcase/math"],
    ["the home page in German", "de", "", "/1.0/de"],
    ["a page the target version lacks", "en", "showcase/widgets", "/1.0"],
    ["a German page the target version lacks", "de", "showcase/math", "/1.0"],
  ] as const)("%s", (_, locale, pagePath, expected) => {
    expect(versionPageUrl("1.0", manifest, locale, pagePath)).toBe(expected);
  });

  test("a target version without a page list opens its root", () => {
    expect(versionPageUrl("1.0", undefined, "en", "showcase/math")).toBe("/1.0");
  });
});
