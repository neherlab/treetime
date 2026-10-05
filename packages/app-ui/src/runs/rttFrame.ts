import { type AxisFrame, clamp, niceAxis } from "./palette";
import type { PointPayload, RttPoint } from "./RootToTipPlot";

export function plotFrame(points: readonly RttPoint[], fitToModel: boolean): PlotFrame {
  const inModel = points.filter((point) => !point.excluded);
  const basis = fitToModel && inModel.length > 0 ? inModel : points;
  const dates = basis.map((point) => point.date);
  const divs = basis.map((point) => point.div);

  return { x: niceAxis(dates), y: niceAxis([0, ...divs]) };
}

export function placePoints(points: readonly RttPoint[], frame: PlotFrame): PlacedPoint[] {
  return points.map((point) => {
    const x = clamp(point.date, ...frame.x.domain);
    const y = clamp(point.div, ...frame.y.domain);

    return { ...point, x, y, offAxes: x !== point.date || y !== point.div };
  });
}

export interface PlotFrame {
  x: AxisFrame;
  y: AxisFrame;
}

export interface PlacedPoint extends PointPayload {
  x: number;
  y: number;
}
