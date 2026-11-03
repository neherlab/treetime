import { htmlToSvg } from "auspice/src/components/download/domToSvg";

import { registerDom } from "../../__tests__/dom";

const CHAR_WIDTH = 6;

const LINE_HEIGHT = 12;

beforeAll(() => {
  registerDom();
});

afterEach(() => {
  Reflect.deleteProperty(document, "createRange");
  document.body.replaceChildren();
});

describe("sVG text of an HTML caption (one line on a grid of 6 px characters)", () => {
  test.each([
    {
      caption: "a line of plain text",
      parts: [{ text: "Last updated 2026-11-03" }],
      texts: [{ x: 0, text: "Last updated 2026-11-03" }],
    },
    {
      caption: "a map attribution with text after its link",
      parts: [
        { text: "© " },
        { tag: "a", text: "Natural Earth", href: "https://www.naturalearthdata.com" },
        { text: " contributors" },
      ],
      texts: [
        { x: 0, text: "© " },
        { x: 2 * CHAR_WIDTH, text: "Natural Earth" },
        { x: 16 * CHAR_WIDTH, text: "contributors" },
      ],
    },
    {
      caption: "a line with text after an italic element",
      parts: [{ tag: "i", text: "Nature" }, { text: " 2026" }],
      texts: [
        { x: 0, text: "Nature" },
        { x: 7 * CHAR_WIDTH, text: "2026" },
      ],
    },
  ] satisfies readonly { caption: string; parts: readonly CaptionPart[]; texts: readonly SvgText[] }[])(
    "$caption starts each text at its first visible character (SVG text drops a leading space)",
    ({ parts, texts }) => {
      const caption = captionLine(parts);

      expect(svgTexts(htmlToSvg(caption).markup)).toStrictEqual(texts);
    },
  );
});

function captionLine(parts: readonly CaptionPart[]): HTMLElement {
  const caption = document.createElement("div");
  caption.style.fontSize = "10px";
  caption.append(...parts.map(captionNode));
  document.body.append(caption);

  const textStarts = new Map<Node, number>();
  let column = 0;

  const walker = document.createTreeWalker(caption, NodeFilter.SHOW_TEXT);

  for (let node = walker.nextNode(); node !== null; node = walker.nextNode()) {
    textStarts.set(node, column);
    column += node.nodeValue?.length ?? 0;
  }

  Object.assign(caption, { getBoundingClientRect: () => gridRect(0, column) });
  Object.assign(document, { createRange: () => new GridRange(textStarts) });

  return caption;
}

function captionNode(part: CaptionPart): Node {
  if (part.tag === undefined) {
    return document.createTextNode(part.text);
  }

  const element = document.createElement(part.tag);
  element.textContent = part.text;

  if (part.href !== undefined) {
    element.setAttribute("href", part.href);
  }

  return element;
}

interface CaptionPart {
  tag?: "a" | "i";
  text: string;
  href?: string;
}

function svgTexts(markup: string): SvgText[] {
  return [...markup.matchAll(/<text x="(?<x>[^"]+)"[^>]*>(?<text>[^<]*)<\/text>/gu)].map((match) => ({
    x: Number(match.groups?.["x"]),
    text: (match.groups?.["text"] ?? "").replaceAll("&amp;", "&"),
  }));
}

interface SvgText {
  x: number;
  text: string;
}

class GridRange {
  readonly #textStarts: ReadonlyMap<Node, number>;
  #node: Node | undefined;
  #start = 0;
  #end = 0;

  constructor(textStarts: ReadonlyMap<Node, number>) {
    this.#textStarts = textStarts;
  }

  selectNodeContents(node: Node): void {
    this.#node = node;
    this.#start = 0;
    this.#end = node.nodeValue?.length ?? 0;
  }

  setStart(node: Node, offset: number): void {
    this.#node = node;
    this.#start = offset;
  }

  setEnd(node: Node, offset: number): void {
    this.#node = node;
    this.#end = offset;
  }

  getClientRects(): GridRect[] {
    const column = this.#node === undefined ? undefined : this.#textStarts.get(this.#node);

    return column === undefined || this.#end <= this.#start
      ? []
      : [gridRect(column + this.#start, this.#end - this.#start)];
  }
}

function gridRect(column: number, length: number): GridRect {
  const left = column * CHAR_WIDTH;
  const width = length * CHAR_WIDTH;

  return { left, top: 0, right: left + width, bottom: LINE_HEIGHT, width, height: LINE_HEIGHT, x: left, y: 0 };
}

interface GridRect {
  left: number;
  top: number;
  right: number;
  bottom: number;
  width: number;
  height: number;
  x: number;
  y: number;
}
