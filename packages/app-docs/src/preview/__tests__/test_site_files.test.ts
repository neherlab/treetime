import { siteFile } from "../site-files";

const FOLDERS = { base: "/dev", build: "/site/dist", root: "/site/dist-root", locales: ["en", "de"] };

describe("site files of the local docs server", () => {
  test.each([
    ["the site root page", "/", ["/site/dist-root/index.html"]],
    ["the version list", "/versions.json", ["/site/dist-root/versions.json"]],
    ["the page that leads to the English home page of the version", "/dev", ["/site/dist/index.html"]],
    ["the English home page", "/dev/en", ["/site/dist/en", "/site/dist/en.html"]],
    [
      "a page without its extension",
      "/dev/en/showcase/math",
      ["/site/dist/en/showcase/math", "/site/dist/en/showcase/math.html"],
    ],
    ["the German home page", "/dev/de", ["/site/dist/de", "/site/dist/de.html"]],
    [
      "a hashed asset",
      "/dev/assets/a-12345678.js",
      ["/site/dist/assets/a-12345678.js", "/site/dist/assets/a-12345678.js.html"],
    ],
  ])("serves %s like the docs host", (_, path, candidates) => {
    const target = siteFile(path, FOLDERS);

    expect(target.kind === "file" ? target.candidates : []).toStrictEqual(candidates);
  });

  test.each([
    ["an English page", "/dev/en/showcase/missing", "/site/dist/en/404.html"],
    ["a German page", "/dev/de/showcase/missing", "/site/dist/de/404.html"],
    ["a path outside every language", "/dev/showcase/missing", "/site/dist/404.html"],
  ])("answers a missing %s with the not-found page of its locale", (_, path, page) => {
    const target = siteFile(path, FOLDERS);

    expect(target.kind === "file" ? target.notFound : undefined).toBe(page);
  });

  test("never serves a path that climbs out of the site", () => {
    expect(siteFile("/dev/../../etc/passwd", FOLDERS)).toStrictEqual({ kind: "missing", page: "/site/dist/404.html" });
  });
});
