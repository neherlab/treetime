import { redirectPage } from "../../lib/redirect-page";
import { manifestPages, pageUrl } from "../built-pages";

describe("built pages", () => {
  test.each([
    ["en.html", "/dev/en"],
    ["de.html", "/dev/de"],
    ["en/showcase/math.html", "/dev/en/showcase/math"],
    ["de/showcase/math.html", "/dev/de/showcase/math"],
  ])("the search index lists %s under %s", (file, expected) => {
    expect(pageUrl("/dev", file)).toBe(expected);
  });

  test("the page manifest lists the pages of each locale relative to its root, without the not-found and redirect pages", () => {
    const files = [
      "404.html",
      "de.html",
      "de/404.html",
      "de/showcase/math.html",
      "en.html",
      "en/404.html",
      "en/showcase/math.html",
      "index.html",
    ];

    expect(manifestPages(files)).toStrictEqual({ en: ["", "showcase/math"], de: ["", "showcase/math"] });
  });

  test("the site root page sends the reader to the English home page of the preferred version", () => {
    const html = redirectPage("https://docs.treetime.neherlab.click", "/dev/en");

    expect([
      html.includes('<meta http-equiv="refresh" content="0; url=/dev/en">'),
      html.includes('<link rel="canonical" href="https://docs.treetime.neherlab.click/dev/en">'),
      html.includes('<a href="/dev/en">'),
      html.includes(`<meta http-equiv="content-security-policy" content="default-src 'none'">`),
    ]).toStrictEqual([true, true, true, true]);
  });
});
