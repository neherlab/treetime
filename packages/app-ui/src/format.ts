import { DateTime } from "luxon";

const SIGNIFICANT_DIGITS = 3;

const MILLISECONDS_PER_DAY = 86_400_000;

export function formatBytes(bytes: number): string {
  if (bytes < 1024) {
    return `${bytes} B`;
  }

  if (bytes < 1024 * 1024) {
    return `${(bytes / 1024).toFixed(1)} kB`;
  }

  return `${(bytes / 1024 / 1024).toFixed(1)} MB`;
}

export function formatDuration(seconds: number): string {
  if (seconds < 1) {
    return `${Math.round(seconds * 1000)} ms`;
  }

  if (seconds < 60) {
    return `${seconds.toFixed(1)} s`;
  }

  return `${Math.floor(seconds / 60)} min ${Math.round(seconds % 60)} s`;
}

export function dayLabel(timestamp: string, now: DateTime): string {
  const day = DateTime.fromISO(timestamp).setZone(now.zone);

  if (!day.isValid) {
    return "Unknown date";
  }

  if (day.hasSame(now, "day")) {
    return "Today";
  }

  if (day.hasSame(now.minus({ days: 1 }), "day")) {
    return "Yesterday";
  }

  return day.toFormat("d LLL yyyy");
}

export function headlineText(headline: Readonly<Record<string, number | string>>): string {
  const rootDate = headline["root_date"];

  if (rootDate !== undefined) {
    return Number.isFinite(Number(rootDate)) ? formatMonthYear(Number(rootDate)) : `root ${String(rootDate)}`;
  }

  const rate = headline["clock_rate"];

  if (rate !== undefined) {
    return Number.isFinite(Number(rate)) ? `rate ${formatRate(Number(rate))}` : `rate ${String(rate)}`;
  }

  return "";
}

function formatMonthYear(year: number): string {
  return decimalYearToDate(year).toFormat("LLL yyyy");
}

export function formatDecimalDate(year: number): string {
  return decimalYearToDate(year).toFormat("yyyy-MM-dd");
}

export function daysBetween(from: number, to: number): number {
  return decimalYearToDate(to).diff(decimalYearToDate(from), "days").days;
}

export function formatRate(rate: number): string {
  return rate.toExponential(SIGNIFICANT_DIGITS - 1);
}

export function formatSignedDays(days: number): string {
  const rounded = Math.round(days);

  return `${rounded > 0 ? "+" : ""}${rounded} d`;
}

function decimalYearToDate(year: number): DateTime {
  const whole = Math.floor(year);
  const start = DateTime.utc(whole, 1, 1);
  const days = start.daysInYear * (year - whole);

  return start.plus({ milliseconds: Math.round(days * MILLISECONDS_PER_DAY) });
}
