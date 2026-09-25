import { DateTime } from "luxon";
import { describe, expect, test } from "vitest";

import { dayLabel, formatBytes, formatDuration, headlineText } from "../format";

describe("format", () => {
  test("a decimal root date maps to its calendar month", () => {
    expect([headlineText({ root_date: 2015.5 }), headlineText({ root_date: 2016 })]).toStrictEqual([
      "Jul 2015",
      "Jan 2016",
    ]);
  });

  test("sizes use binary kilobytes and megabytes", () => {
    expect([formatBytes(512), formatBytes(2048), formatBytes(3 * 1024 * 1024)]).toStrictEqual([
      "512 B",
      "2.0 kB",
      "3.0 MB",
    ]);
  });

  test("durations switch units at one second and one minute", () => {
    expect([formatDuration(0.25), formatDuration(12.34), formatDuration(125)]).toStrictEqual([
      "250 ms",
      "12.3 s",
      "2 min 5 s",
    ]);
  });

  test("run days are today, yesterday or a date", () => {
    const now = DateTime.fromISO("2026-09-25T12:00:00Z", { zone: "utc" });

    expect([
      dayLabel("2026-09-25T08:00:00Z", now),
      dayLabel("2026-09-24T08:00:00Z", now),
      dayLabel("2026-09-01T08:00:00Z", now),
    ]).toStrictEqual(["Today", "Yesterday", "1 Sep 2026"]);
  });

  test("the headline shows the root date of a time tree before the rate", () => {
    expect(headlineText({ root_date: 2013.9, clock_rate: 0.0008 })).toStrictEqual("Nov 2013");
  });

  test("the headline shows the clock rate of a clock run", () => {
    expect(headlineText({ clock_rate: 0.00081234, r_squared: 0.9 })).toStrictEqual("rate 8.12e-4");
  });

  test("a non-finite headline value is shown as written", () => {
    expect(headlineText({ clock_rate: "nan" })).toStrictEqual("rate nan");
  });
});
