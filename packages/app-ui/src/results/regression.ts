export interface Line {
  slope: number;
  intercept: number;
}

export interface Point {
  x: number;
  y: number;
}

export function leastSquares(points: readonly Point[]): Line | undefined {
  if (points.length < 2) {
    return undefined;
  }

  const meanX = mean(points.map((point) => point.x));
  const meanY = mean(points.map((point) => point.y));
  const sxx = points.reduce((sum, point) => sum + (point.x - meanX) ** 2, 0);
  const sxy = points.reduce((sum, point) => sum + (point.x - meanX) * (point.y - meanY), 0);

  if (sxx === 0) {
    return undefined;
  }

  const slope = sxy / sxx;

  return { slope, intercept: meanY - slope * meanX };
}

function mean(values: readonly number[]): number {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}
