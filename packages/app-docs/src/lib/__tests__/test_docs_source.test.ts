import { createDocsSource, type DocsPageData, isFallbackPage } from "../docs-source";

describe("docs source", () => {
  test("a German path without a German file gives the English page, marked as a fallback", () => {
    const source = createDocsSource(FILES, "/dev");
    const math = source.getPage(["showcase", "math"], "de");

    expect([math?.data.title, math?.url, math === undefined ? undefined : isFallbackPage(math, "de")]).toStrictEqual([
      "Math",
      "/dev/de/showcase/math",
      true,
    ]);
  });

  test("a German path with a German file gives the German page", () => {
    const source = createDocsSource(FILES, "/dev");
    const home = source.getPage([], "de");

    expect([home?.data.title, home?.url, home === undefined ? undefined : isFallbackPage(home, "de")]).toStrictEqual([
      "Startseite",
      "/dev/de",
      false,
    ]);
  });

  test("an English page lies below the English language segment", () => {
    const source = createDocsSource(FILES, "/dev");
    const math = source.getPage(["showcase", "math"], "en");

    expect([math?.url, math === undefined ? undefined : isFallbackPage(math, "en")]).toStrictEqual([
      "/dev/en/showcase/math",
      false,
    ]);
  });

  test("the page tree links to the docs urls of its locale, in the order of meta.json", () => {
    const source = createDocsSource(FILES, "/dev");
    const folder = source.getPageTree("de").children.find((node) => node.type === "folder");

    expect(
      folder?.type === "folder" ? folder.children.map((node) => (node.type === "page" ? node.url : node.type)) : [],
    ).toStrictEqual(["/dev/de/showcase/typography", "/dev/de/showcase/math"]);
  });
});

const FILES = {
  pages: [
    { path: "en/index.mdx", data: pageData("Home", "en/index.mdx") },
    { path: "de/index.mdx", data: pageData("Startseite", "de/index.mdx") },
    { path: "en/showcase/math.mdx", data: pageData("Math", "en/showcase/math.mdx") },
    { path: "en/showcase/typography.mdx", data: pageData("Typography", "en/showcase/typography.mdx") },
  ],
  metas: [
    { path: "en/meta.json", data: { pages: ["index", "showcase"] } },
    { path: "de/meta.json", data: { pages: ["index", "showcase"] } },
    { path: "en/showcase/meta.json", data: { title: "Showcase", pages: ["typography", "math"] } },
  ],
};

function pageData(title: string, file: string): DocsPageData {
  return {
    title,
    template: "doc",
    tableOfContents: true,
    body: { kind: "entry", key: `docs:${file}` },
    editPath: file,
  };
}
