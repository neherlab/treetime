export function redirectPage(site: string, target: string): string {
  const canonical = new URL(target, site).href;

  return [
    "<!doctype html>",
    '<html lang="en">',
    "<head>",
    '<meta charset="utf-8">',
    `<meta http-equiv="content-security-policy" content="default-src 'none'">`,
    "<title>TreeTime docs</title>",
    `<meta http-equiv="refresh" content="0; url=${target}">`,
    `<link rel="canonical" href="${canonical}">`,
    "</head>",
    "<body>",
    `<p><a href="${target}">TreeTime docs</a></p>`,
    "</body>",
    "</html>",
    "",
  ].join("\n");
}
