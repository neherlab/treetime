const WORK_AREA_SHARE = 0.9;

const MAX_FIRST_SIZE = { width: 1600, height: 1000 } as const;

export interface Rectangle {
  x: number;
  y: number;
  width: number;
  height: number;
}

export interface Size {
  width: number;
  height: number;
}

export function firstWindowBounds(workArea: Rectangle, minSize: Size): Rectangle {
  const width = firstLength(workArea.width, MAX_FIRST_SIZE.width, minSize.width);
  const height = firstLength(workArea.height, MAX_FIRST_SIZE.height, minSize.height);

  return {
    x: workArea.x + Math.floor((workArea.width - width) / 2),
    y: workArea.y + Math.floor((workArea.height - height) / 2),
    width,
    height,
  };
}

function firstLength(available: number, max: number, min: number): number {
  return Math.round(Math.min(Math.max(Math.min(available * WORK_AREA_SHARE, max), min), available));
}
