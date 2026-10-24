import { execFileSync } from "node:child_process";
import { createHash } from "node:crypto";
import { mkdir, mkdtemp, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

import * as z from "zod";

const NATURAL_EARTH = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/v5.1.2/geojson";

const FONTS = "https://github.com/openmaptiles/fonts";

const DOWNLOADS = {
  countries: {
    url: `${NATURAL_EARTH}/ne_50m_admin_0_countries.geojson`,
    sha256: "3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb",
  },
  borders: {
    url: `${NATURAL_EARTH}/ne_50m_admin_0_boundary_lines_land.geojson`,
    sha256: "2faac4f6b34386f3d21b6e018cf151f241f00e5c936d44dd17d7d9bfb147fa48",
  },
  fonts: {
    url: `${FONTS}/releases/download/v2.0/noto-sans.zip`,
    sha256: "d117316544b43a5dde7ee761b36e17701e9f85574e181d76a74814240fdbaf34",
  },
  fontLicense: {
    url: "https://raw.githubusercontent.com/openmaptiles/fonts/v2.0/noto-sans/LICENSE",
    sha256: "6a73f9541c2de74158c0e7cf6b0a58ef774f5a780bf191f2d7ec9cc53efe2bf2",
  },
} as const;

const DOWNLOAD_TIMEOUT_MS = 300_000;

const GLYPHS = "Noto Sans Bold/0-255.pbf";

const BASEMAP_DIR = fileURLToPath(new URL("../public/basemap", import.meta.url));

const COORDINATE_STEP = 100;

const UNZIP_MAX_BYTES = 10_000_000;

const position = z.array(z.number());

const line = z.array(position);

const geometrySchema = z.discriminatedUnion("type", [
  z.object({ type: z.literal("LineString"), coordinates: line }),
  z.object({ type: z.literal("MultiLineString"), coordinates: z.array(line) }),
  z.object({ type: z.literal("Polygon"), coordinates: z.array(line) }),
  z.object({ type: z.literal("MultiPolygon"), coordinates: z.array(z.array(line)) }),
]);

const countriesSchema = z.object({
  type: z.literal("FeatureCollection"),
  features: z.array(
    z.object({
      properties: z.object({ NAME: z.string(), LABEL_X: z.number(), LABEL_Y: z.number(), LABELRANK: z.number() }),
      geometry: geometrySchema,
    }),
  ),
});

const bordersSchema = z.object({
  type: z.literal("FeatureCollection"),
  features: z.array(z.object({ geometry: geometrySchema })),
});

if (import.meta.main) {
  await writeBasemap();
}

export function countryOutlines(countries: Countries): FeatureCollection {
  return featureCollection(
    countries.features.map(({ properties, geometry }) => feature({ name: properties.NAME }, roundGeometry(geometry))),
  );
}

export function landBorders(borders: Borders): FeatureCollection {
  return featureCollection(borders.features.map(({ geometry }) => feature({}, roundGeometry(geometry))));
}

export function countryLabels(countries: Countries): FeatureCollection {
  return featureCollection(
    countries.features.map(({ properties }) =>
      feature(
        { name: properties.NAME, rank: labelRank(properties.LABELRANK) },
        { type: "Point", coordinates: roundPosition([properties.LABEL_X, properties.LABEL_Y]) },
      ),
    ),
  );
}

export function labelRank(naturalEarthRank: number): 1 | 2 | 3 {
  if (naturalEarthRank <= 2) {
    return 1;
  }

  return naturalEarthRank === 3 ? 2 : 3;
}

export function parseCountries(text: string): Countries {
  return countriesSchema.parse(JSON.parse(text));
}

export function parseBorders(text: string): Borders {
  return bordersSchema.parse(JSON.parse(text));
}

export type Countries = z.infer<typeof countriesSchema>;

export type Borders = z.infer<typeof bordersSchema>;

export interface FeatureCollection {
  type: "FeatureCollection";
  features: Feature[];
}

interface Feature {
  type: "Feature";
  properties: Readonly<Record<string, string | number>>;
  geometry: Geometry | Point;
}

type Geometry = z.infer<typeof geometrySchema>;

interface Point {
  type: "Point";
  coordinates: number[];
}

async function writeBasemap(): Promise<void> {
  const countries = parseCountries(await downloadText(DOWNLOADS.countries));
  const borders = parseBorders(await downloadText(DOWNLOADS.borders));
  const glyphs = await zipEntry(await download(DOWNLOADS.fonts), GLYPHS);
  const license = await downloadText(DOWNLOADS.fontLicense);

  const glyphsPath = join(BASEMAP_DIR, "fonts", GLYPHS);

  await mkdir(dirname(glyphsPath), { recursive: true });
  await writeFile(join(BASEMAP_DIR, "countries.json"), JSON.stringify(countryOutlines(countries)));
  await writeFile(join(BASEMAP_DIR, "borders.json"), JSON.stringify(landBorders(borders)));
  await writeFile(join(BASEMAP_DIR, "labels.json"), JSON.stringify(countryLabels(countries)));
  await writeFile(glyphsPath, glyphs);
  await writeFile(join(BASEMAP_DIR, "fonts", "LICENSE.txt"), license);
}

async function downloadText(source: Download): Promise<string> {
  return new TextDecoder().decode(await download(source));
}

async function download({ url, sha256 }: Download): Promise<Uint8Array> {
  const response = await fetch(url, { signal: AbortSignal.timeout(DOWNLOAD_TIMEOUT_MS) });

  if (!response.ok) {
    throw new Error(`Downloading ${url} failed: HTTP ${response.status}`);
  }

  const bytes = new Uint8Array(await response.arrayBuffer());
  const digest = createHash("sha256").update(bytes).digest("hex");

  if (digest !== sha256) {
    throw new Error(`${url} has SHA-256 ${digest}, expected ${sha256}`);
  }

  return bytes;
}

interface Download {
  url: string;
  sha256: string;
}

async function zipEntry(archive: Uint8Array, entry: string): Promise<Uint8Array> {
  const dir = await mkdtemp(join(tmpdir(), "treetime-basemap-"));

  try {
    const path = join(dir, "archive.zip");

    await writeFile(path, archive);

    return new Uint8Array(execFileSync("unzip", ["-p", path, entry], { maxBuffer: UNZIP_MAX_BYTES }));
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
}

function featureCollection(features: Feature[]): FeatureCollection {
  return { type: "FeatureCollection", features };
}

function feature(properties: Feature["properties"], geometry: Feature["geometry"]): Feature {
  return { type: "Feature", properties, geometry };
}

function roundGeometry(geometry: Geometry): Geometry {
  switch (geometry.type) {
    case "LineString":
      return { type: geometry.type, coordinates: roundLine(geometry.coordinates) };
    case "MultiLineString":
    case "Polygon":
      return { type: geometry.type, coordinates: geometry.coordinates.map(roundLine) };
    case "MultiPolygon":
      return { type: geometry.type, coordinates: geometry.coordinates.map((polygon) => polygon.map(roundLine)) };
    default:
      return assertNever(geometry);
  }
}

function roundLine(positions: number[][]): number[][] {
  return positions.map(roundPosition);
}

function roundPosition(degrees: number[]): number[] {
  return degrees.map((value) => Math.round(value * COORDINATE_STEP) / COORDINATE_STEP);
}

function assertNever(value: never): never {
  throw new Error(`Unhandled geometry: ${JSON.stringify(value)}`);
}
