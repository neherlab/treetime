import type { WebContents } from "electron";

export function contentSecurityPolicy(devServer: boolean): string {
  const scriptSources = devServer ? "'self' 'unsafe-inline'" : "'self'";

  return [
    "default-src 'self'",
    `script-src ${scriptSources}`,
    "style-src 'self' 'unsafe-inline'",
    "img-src 'self' data: blob:",
    "font-src 'self' data:",
    "connect-src 'self'",
    "worker-src 'self' blob:",
    "object-src 'none'",
    "base-uri 'none'",
    "form-action 'none'",
  ].join("; ");
}

export function isAppUrl(url: string, appUrl: string): boolean {
  const candidate = parseUrl(url);
  const app = parseUrl(appUrl);

  if (candidate === undefined || app === undefined || candidate.protocol !== app.protocol) {
    return false;
  }

  if (app.protocol === "file:") {
    return candidate.pathname === app.pathname;
  }

  return candidate.origin === app.origin;
}

export interface SenderFrame {
  url: string;
  parent: unknown;
}

export function isTrustedSender(frame: SenderFrame | null, appUrl: string): boolean {
  return frame !== null && frame.parent === null && isAppUrl(frame.url, appUrl);
}

export function isExternalLink(url: string): boolean {
  const parsed = parseUrl(url);

  return parsed !== undefined && (parsed.protocol === "https:" || parsed.protocol === "http:");
}

export function confineNavigation(contents: WebContents, appUrl: string, openExternal: (url: string) => void): void {
  contents.setWindowOpenHandler(({ url }) => {
    if (isExternalLink(url)) {
      openExternal(url);
    }

    return { action: "deny" };
  });

  contents.on("will-navigate", (event, url) => {
    if (!isAppUrl(url, appUrl)) {
      event.preventDefault();
    }
  });

  contents.on("will-attach-webview", (event) => {
    event.preventDefault();
  });
}

function parseUrl(url: string): URL | undefined {
  try {
    return new URL(url);
  } catch {
    return undefined;
  }
}
