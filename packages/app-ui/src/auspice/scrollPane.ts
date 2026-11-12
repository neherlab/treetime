const SCROLLING = new Set(["auto", "scroll"]);

export function scrollAncestor(element: Element | null): HTMLElement | undefined {
  for (let node = element?.parentElement ?? null; node !== null; node = node.parentElement) {
    if (SCROLLING.has(getComputedStyle(node).overflowY)) {
      return node;
    }
  }

  return undefined;
}

export function observeHeight(element: HTMLElement | undefined, onHeight: (height: number) => void): () => void {
  const measure = () => onHeight(element?.clientHeight ?? window.innerHeight);

  measure();

  if (element === undefined) {
    window.addEventListener("resize", measure);

    return () => window.removeEventListener("resize", measure);
  }

  const observer = new ResizeObserver(measure);
  observer.observe(element);

  return () => observer.disconnect();
}
