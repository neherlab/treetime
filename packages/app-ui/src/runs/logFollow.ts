export const FOLLOW_THRESHOLD_PX = 8;

export function isScrolledToEnd({
  scrollTop,
  clientHeight,
  scrollHeight,
}: {
  scrollTop: number;
  clientHeight: number;
  scrollHeight: number;
}): boolean {
  return scrollTop + clientHeight >= scrollHeight - FOLLOW_THRESHOLD_PX;
}
