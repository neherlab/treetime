function isAppUrl(url: string, appUrl: string): boolean {
  const candidate = parseUrl(url);
  const app = parseUrl(appUrl);

  return (
    candidate !== undefined &&
    app !== undefined &&
    candidate.protocol === app.protocol &&
    candidate.host === app.host &&
    candidate.host !== ""
  );
}

interface SenderFrame {
  url: string;
  parent: unknown;
}

export function isTrustedSender(frame: SenderFrame | null, appUrl: string): boolean {
  return frame !== null && frame.parent === null && isAppUrl(frame.url, appUrl);
}

function isExternalLink(url: string): boolean {
  const parsed = parseUrl(url);

  return parsed !== undefined && (parsed.protocol === "https:" || parsed.protocol === "http:");
}

export function confineNavigation(
  contents: NavigableContents,
  appUrl: string,
  openExternal: (url: string) => void,
): void {
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

interface PreventableEvent {
  preventDefault(): void;
}

interface NavigableContents {
  setWindowOpenHandler(handler: (details: { url: string }) => { action: "deny" }): void;
  on(event: "will-navigate", listener: (event: PreventableEvent, url: string) => void): void;
  on(event: "will-attach-webview", listener: (event: PreventableEvent) => void): void;
}

function parseUrl(url: string): URL | undefined {
  try {
    return new URL(url);
  } catch {
    return undefined;
  }
}
