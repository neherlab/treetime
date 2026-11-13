export function onPageHide(target: EventTarget, listener: () => void): () => void {
  target.addEventListener("pagehide", listener, { once: true });

  return () => target.removeEventListener("pagehide", listener);
}
