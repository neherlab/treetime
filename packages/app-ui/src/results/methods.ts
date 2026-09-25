import { formatDecimalDate, formatRate } from "../format";
import { isNumber, type JsonObject } from "../settings/json";
import type { TimetreeEstimates } from "./estimates";

export const CITATION =
  "Sagulenko P, Puller V, Neher RA. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4 (2018), vex042.";

export const CITATION_DOI = "https://doi.org/10.1093/ve/vex042";

export function coalescentPrior(config: JsonObject): string {
  const tc = config["coalescent"];

  if (config["coalescent_skyline"] === true) {
    const points = config["skyline_n_points"];
    const stiffness = config["skyline_stiffness"];

    const details = [isNumber(points) ? `${points} points` : "", isNumber(stiffness) ? `stiffness ${stiffness}` : ""];

    return ["Skyline", ...details.filter((detail) => detail !== "")].join(", ");
  }

  if (config["coalescent_opt"] === true) {
    return "Constant size, optimized Tc";
  }

  if (isNumber(tc)) {
    return `Constant size, Tc = ${tc} years`;
  }

  return "None";
}

export function relaxedClock(config: JsonObject): readonly [number, number] | undefined {
  const relax = config["relax"];

  return Array.isArray(relax) && relax.length === 2 && isNumber(relax[0]) && isNumber(relax[1])
    ? [relax[0], relax[1]]
    : undefined;
}

export function timetreeMethods(version: string, config: JsonObject, estimates: TimetreeEstimates): string {
  const sentences = [
    `A time-scaled phylogeny of ${estimates.samples} samples was inferred with TreeTime ${version} (timetree command).`,
    rateSentence(config, estimates),
    filterSentence(config, estimates),
    priorSentence(config),
    relaxSentence(config),
    rootSentence(estimates),
    `Please cite: ${CITATION}`,
  ];

  return sentences.filter((sentence) => sentence !== "").join(" ");
}

function rateSentence(config: JsonObject, estimates: TimetreeEstimates): string {
  if (estimates.rate === undefined) {
    return "";
  }

  if (estimates.rateFixed) {
    const std = config["clock_std_dev"];

    const spread = isNumber(std) ? ` with standard deviation ${formatRate(std)}` : "";

    return `The clock rate was fixed at ${formatRate(estimates.rate)} substitutions per site per year${spread}.`;
  }

  const spread = estimates.rateStd === undefined ? "" : ` (standard deviation ${formatRate(estimates.rateStd)})`;

  return `The clock rate was estimated at ${formatRate(estimates.rate)} substitutions per site per year${spread}.`;
}

function filterSentence(config: JsonObject, estimates: TimetreeEstimates): string {
  const threshold = config["clock_filter"];

  const filter =
    isNumber(threshold) && threshold > 0
      ? `Samples whose root-to-tip residual exceeded ${threshold} interquartile distances were treated as clock outliers. `
      : "";

  return `${filter}${estimates.excludedSamples} of ${estimates.samples} samples had no usable date or were clock outliers and did not constrain the clock model.`;
}

function priorSentence(config: JsonObject): string {
  const prior = coalescentPrior(config);

  return prior === "None" ? "" : `A coalescent prior was used (${prior.toLowerCase()}).`;
}

function relaxSentence(config: JsonObject): string {
  const relax = relaxedClock(config);

  return relax === undefined ? "" : `A relaxed clock was used (slack ${relax[0]}, coupling ${relax[1]}).`;
}

function rootSentence(estimates: TimetreeEstimates): string {
  if (estimates.rootDate === undefined) {
    return "";
  }

  const interval = estimates.rootInterval;

  const range =
    interval === undefined
      ? ""
      : ` (90% interval ${formatDecimalDate(interval[0])} to ${formatDecimalDate(interval[1])})`;

  return `The root was dated to ${formatDecimalDate(estimates.rootDate)}${range}.`;
}
