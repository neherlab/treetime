import { pageView } from "../page-view";

describe("page view", () => {
  test("the table of contents lists headings of levels 2 and 3 without the hidden footnote label", () => {
    const view = pageView({
      page: {
        url: "/dev/showcase/typography",
        path: "showcase/typography.mdx",
        slugs: ["showcase", "typography"],
        data: {
          title: "Typography",
          template: "doc",
          tableOfContents: true,
          entryId: "showcase/typography",
          file: "src/content/docs/showcase/typography.mdx",
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
});
