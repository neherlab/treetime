import { createHash } from "node:crypto";

import { Window } from "happy-dom";

const JAVASCRIPT_TYPES = new Set([
  "",
  "module",
  "text/javascript",
  "application/javascript",
  "application/ecmascript",
  "text/ecmascript",
]);

export function cspViolations(page: string, html: string): string[] {
  const window = new Window();

  try {
    const document = new window.DOMParser().parseFromString(html, "text/html");

    const policy = document.querySelector('meta[http-equiv="content-security-policy" i]')?.getAttribute("content");

    if (policy === null || policy === undefined) {
      return [`${page}: the page has no content security policy meta element`];
    }

    const allowed = scriptSources(policy);

    return [...document.querySelectorAll("script:not([src])")].flatMap((script) => {
      const executable = JAVASCRIPT_TYPES.has((script.getAttribute("type") ?? "").trim().toLowerCase());
      const hash = `'sha256-${createHash("sha256").update(script.textContent).digest("base64")}'`;

      return executable && !allowed.has(hash)
        ? [`${page}: the inline script with hash ${hash} is missing from script-src`]
        : [];
    });
  } finally {
    void window.happyDOM.close();
  }
}

function scriptSources(policy: string): Set<string> {
  return new Set(
    policy.split(";").flatMap((directive) => {
      const [name, ...sources] = directive.trim().split(/\s+/u);

      return name === "script-src" || name === "script-src-elem" ? sources : [];
    }),
  );
}
