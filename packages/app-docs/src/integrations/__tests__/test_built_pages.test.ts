import { manifestPages, pageUrl } from "../built-pages";
import { rootRedirectPage } from "../site-root";

describe("built pages", () => {
  test.each([
    ["index.html", "/dev"],
    ["de.html", "/dev/de"],
    ["showcase/math.html", "/dev/showcase/math"],
    ["de/showcase/math.html", "/dev/de/showcase/math"],
  ])("the search index lists %s under %s", (file, expected) => {
    expect(pageUrl("/dev", file)).toBe(expected);
  });

  test("the page manifest lists the pages of each locale relative to its root, without the not-found pages", () => {
    const files = ["404.html", "de.html", "de/404.html", "de/showcase/math.html", "index.html", "showcase/math.html"];

    expect(manifestPages(files)).toStrictEqual({ en: ["", "showcase/math"], de: ["", "showcase/math"] });
  });

  test("the site root page sends the reader to the preferred version", () => {
    const html = rootRedirectPage("https://docs.treetime.neherlab.click", "dev");

    expect([
      html.includes('<meta http-equiv="refresh" content="0; url=/dev">'),
      html.includes('<link rel="canonical" href="https://docs.treetime.neherlab.click/dev">'),
      html.includes('<a href="/dev">'),
    ]).toStrictEqual([true, true, true]);
  });
});
