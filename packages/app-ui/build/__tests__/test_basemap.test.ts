import { countryLabels, countryOutlines, labelRank, landBorders, parseBorders, parseCountries } from "../basemap";

const COUNTRIES = parseCountries(
  JSON.stringify({
    type: "FeatureCollection",
    features: [
      {
        type: "Feature",
        properties: { NAME: "Peru", ADMIN: "Peru", LABEL_X: -72.900016, LABEL_Y: -12.976679, LABELRANK: 2 },
        geometry: {
          type: "Polygon",
          coordinates: [
            [
              [-69.896, -4.228],
              [-70.004, -4.3249],
              [-69.9, -4.2312],
              [-69.896, -4.228],
            ],
          ],
        },
      },
    ],
  }),
);

describe("basemap from Natural Earth", () => {
  test("country outlines keep the name and coordinates rounded to 0.01 degree", () => {
    expect(countryOutlines(COUNTRIES)).toStrictEqual({
      type: "FeatureCollection",
      features: [
        {
          type: "Feature",
          properties: { name: "Peru" },
          geometry: {
            type: "Polygon",
            coordinates: [
              [
                [-69.9, -4.23],
                [-70, -4.32],
                [-69.9, -4.23],
                [-69.9, -4.23],
              ],
            ],
          },
        },
      ],
    });
  });

  test("land borders keep no properties", () => {
    const borders = parseBorders(
      JSON.stringify({
        type: "FeatureCollection",
        features: [
          {
            type: "Feature",
            properties: { FEATURECLA: "International boundary (verify)", MIN_ZOOM: 0 },
            geometry: {
              type: "MultiLineString",
              coordinates: [
                [
                  [10.123, 47.456],
                  [10.2, 47.5],
                ],
              ],
            },
          },
        ],
      }),
    );

    expect(landBorders(borders)).toStrictEqual({
      type: "FeatureCollection",
      features: [
        {
          type: "Feature",
          properties: {},
          geometry: {
            type: "MultiLineString",
            coordinates: [
              [
                [10.12, 47.46],
                [10.2, 47.5],
              ],
            ],
          },
        },
      ],
    });
  });

  test("a country label sits at the label point of the country, with its rank", () => {
    expect(countryLabels(COUNTRIES)).toStrictEqual({
      type: "FeatureCollection",
      features: [
        {
          type: "Feature",
          properties: { name: "Peru", rank: 1 },
          geometry: { type: "Point", coordinates: [-72.9, -12.98] },
        },
      ],
    });
  });

  test.each([
    { naturalEarth: 0, rank: 1 },
    { naturalEarth: 2, rank: 1 },
    { naturalEarth: 3, rank: 2 },
    { naturalEarth: 4, rank: 3 },
    { naturalEarth: 10, rank: 3 },
  ])("label rank $naturalEarth of Natural Earth is label rank $rank", ({ naturalEarth, rank }) => {
    expect(labelRank(naturalEarth)).toBe(rank);
  });

  test("a country without a label point is refused", () => {
    expect(() =>
      parseCountries(
        JSON.stringify({
          type: "FeatureCollection",
          features: [{ properties: { NAME: "Peru" }, geometry: { type: "Polygon", coordinates: [] } }],
        }),
      ),
    ).toThrow("LABEL_X");
  });
});
