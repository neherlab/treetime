import type { RunHeadline, zYearDate } from "@neherlab/app-contracts";
import { DateTime } from "luxon";
import type { z } from "zod";

import { fromJsonFloat } from "./results/numbers";

const SIGNIFICANT_DIGITS = 3;

const YEAR_DECIMALS = 2;

const LUXON_YEARS = { first: 0, end: 10_000 };

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

export function headlineText(headline: RunHeadline): string {
  const rootDate = headline.root_date;

  if (rootDate !== undefined) {
    return rootDate.date !== undefined && rootDate.year >= LUXON_YEARS.first && rootDate.year < LUXON_YEARS.end
      ? DateTime.fromISO(rootDate.date, { zone: "utc" }).toFormat("LLL yyyy")
      : yearDateText(rootDate);
  }

  const rate = headline.clock_rate;

  if (rate !== undefined) {
    const value = fromJsonFloat(rate);

    return Number.isFinite(value) ? `rate ${formatRate(value)}` : `rate ${rate}`;
  }

  const recurrent = headline.recurrent_substitutions;

  return recurrent === undefined ? "" : `${recurrent} recurrent`;
}

export function yearDateText(date: z.output<typeof zYearDate>): string {
  return date.date ?? date.year.toFixed(YEAR_DECIMALS);
}

export function formatRate(rate: number): string {
  return rate.toExponential(SIGNIFICANT_DIGITS - 1);
}

export function rSquaredText(rSquared: number | undefined): string | undefined {
  return rSquared === undefined ? undefined : `R² = ${rSquared.toFixed(3)}`;
}

export function formatLevel(level: number): string {
  return `${Math.round(level * 100)}%`;
}

export function formatSignedDays(days: number): string {
  const rounded = Math.round(days);

  return `${rounded > 0 ? "+" : ""}${rounded} d`;
}
