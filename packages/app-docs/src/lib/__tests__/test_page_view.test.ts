import { pageView } from "../page-view";

describe("page view", () => {
  test("the table of contents lists headings of levels 2 and 3 without the hidden footnote label", () => {
    const view = pageView({
      page: {
        url: "/dev/en/showcase/typography",
        path: "en/showcase/typography.mdx",
        slugs: ["showcase", "typography"],
        data: {
          title: "Typography",
          template: "doc",
          tableOfContents: true,
          body: { kind: "entry", key: "docs:en/showcase/typography" },
          editPath: "en/showcase/typography.mdx",
        },
      },
      locale: "en",
      tree: { name: "Docs", children: [] },
      headings: [
        { depth: 2, slug: "lists", text: "Lists" },
        { depth: 3, slug: "nested", text: "Nested" },
        { depth: 4, slug: "deep", text: "Deep" },
        { depth: 2, slug: "footnote-label", text: "Footnotes" },
      ],
      base: "/dev",
    });

    expect(view.tocHeadings.map((heading) => heading.slug)).toStrictEqual(["lists", "nested"]);
  });

  test("the edit link opens the file of the page in the docs folder", () => {
    expect(view({ editPath: "en/showcase/typography.mdx" }).editUrl).toBe(
      "https://github.com/neherlab/treetime/edit/rust/docs/en/showcase/typography.mdx",
    );
  });

  test("a page without a file in the docs folder has no edit link", () => {
    expect(view({ editPath: undefined }).editUrl).toBeUndefined();
  });
});

function view({ editPath }: { editPath: string | undefined }) {
  return pageView({
    page: {
      url: "/dev/en/showcase/typography",
      path: "en/showcase/typography.mdx",
      slugs: ["showcase", "typography"],
      data: {
        title: "Typography",
        template: "doc",
        tableOfContents: true,
        body: { kind: "settings", command: "clock" },
        editPath,
      },
    },
    locale: "en",
    tree: { name: "Docs", children: [] },
    headings: [],
    base: "/dev",
  });
}
