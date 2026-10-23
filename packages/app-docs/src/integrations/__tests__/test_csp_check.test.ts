import { cspViolations } from "../csp-check";

const HASHED_SCRIPT = "document.documentElement.classList.add('dark')";

const OTHER_SCRIPT = "console.log('unlisted')";

describe("content security policy check", () => {
  test("accepts an executable inline script whose hash the policy lists, and a JSON data block", () => {
    expect(
      cspViolations(
        "index.html",
        html([`<script>${HASHED_SCRIPT}</script>`, '<script type="application/json">{"a":1}</script>']),
      ),
    ).toStrictEqual([]);
  });

  test("names the page and the hash of an inline script that the policy lacks", () => {
    expect(cspViolations("index.html", html([`<script type="module">${OTHER_SCRIPT}</script>`]))).toStrictEqual([
      `index.html: the inline script with hash 'sha256-${sha256(OTHER_SCRIPT)}' is missing from script-src`,
    ]);
  });

  test("reports a page without a policy", () => {
    expect(cspViolations("de.html", "<!doctype html><html><head></head><body></body></html>")).toStrictEqual([
      "de.html: the page has no content security policy meta element",
    ]);
  });
});

function html(scripts: string[]): string {
  const policy = `default-src 'self'; script-src 'self' 'sha256-${sha256(HASHED_SCRIPT)}'`;

  return `<!doctype html><html><head><meta http-equiv="content-security-policy" content="${policy}">${scripts.join("")}</head><body></body></html>`;
}

function sha256(text: string): string {
  return new Bun.CryptoHasher("sha256").update(text).digest("base64");
}
