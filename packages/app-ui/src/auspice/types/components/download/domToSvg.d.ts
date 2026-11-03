export declare function htmlToSvg(root: Element, options?: { skip?: (el: Element) => boolean }): SvgBlock;

interface SvgBlock {
  width: number;
  height: number;
  markup: string;
}
