export function rootRedirectPage(site: string, preferred: string): string {
  const target = `/${preferred}`;
  const canonical = new URL(target, site).href;

  return [
    "<!doctype html>",
    '<html lang="en">',
    "<head>",
    '<meta charset="utf-8">',
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
