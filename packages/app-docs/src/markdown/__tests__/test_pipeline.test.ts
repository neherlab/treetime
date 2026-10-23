import { docsMarkdownProcessor } from "../pipeline";

const ENGLISH_PAGE = new URL("file:///docs/src/content/docs/showcase/code.md");

const GERMAN_PAGE = new URL("file:///docs/src/content/docs/de/showcase/code.md");

describe("docs markdown pipeline", () => {
  test("keeps text and leaf directives that no step handles as written", async () => {
    const html = await renderMarkdown("Read `x` as NHX:key=value or ref:main.\n\n::name\n", ENGLISH_PAGE);

    expect(html).toBe("<p>Read <code>x</code> as NHX:key=value or ref:main.</p>\n<p>::name</p>");
  });

  test.each([
    ["an English page", ENGLISH_PAGE, '<p><a href="/dev/showcase/math#eq">math</a></p>'],
    ["a German page", GERMAN_PAGE, '<p><a href="/dev/de/showcase/math#eq">math</a></p>'],
  ])("puts the version base before a root-relative link on %s", async (_, file, expected) => {
    expect(await renderMarkdown("[math](/showcase/math#eq)", file)).toBe(expected);
  });

  test("puts the version base before a file of the public folder without a locale", async () => {
    expect(await renderMarkdown("![logo](/logo.svg)", GERMAN_PAGE)).toBe('<p><img src="/dev/logo.svg" alt="logo"></p>');
  });

  test("stops at a callout type it does not know", async () => {
    await expect(renderMarkdown(":::nte\ntext\n:::\n", ENGLISH_PAGE)).rejects.toThrow(
      "The block `:::nte` is no callout type. Use note, tip, warning, or danger.",
    );
  });
});

async function renderMarkdown(markdown: string, fileURL: URL): Promise<string> {
  const renderer = await docsMarkdownProcessor("/dev").createRenderer({ syntaxHighlight: false });

  return (await renderer.render(markdown, { fileURL })).code;
}
