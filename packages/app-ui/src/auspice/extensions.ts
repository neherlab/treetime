const ATTRIBUTION =
  '<a href="https://www.naturalearthdata.com/" target="_blank" rel="noopener">Natural Earth</a> 5.1.2';

const LABEL_COLOR = "#6b6b6b";

const extensions = {
  get mapTiles(): MapTiles {
    return { style: basemapStyle(), attribution: ATTRIBUTION };
  },
};

export default extensions;

interface MapTiles {
  style: ReturnType<typeof basemapStyle>;
  attribution: string;
}

function basemapStyle() {
  const base = new URL(`${import.meta.env.BASE_URL}basemap/`, document.baseURI).href;

  return {
    version: 8,
    glyphs: `${base}fonts/{fontstack}/{range}.pbf`,
    sources: {
      countries: { type: "geojson", data: `${base}countries.json`, attribution: ATTRIBUTION },
      borders: { type: "geojson", data: `${base}borders.json` },
      labels: { type: "geojson", data: `${base}labels.json` },
    },
    layers: [
      { id: "water", type: "background", paint: { "background-color": "hsl(185, 9%, 81%)" } },
      {
        id: "land",
        type: "fill",
        source: "countries",
        paint: { "fill-color": "hsl(55, 11%, 96%)", "fill-antialias": true },
      },
      {
        id: "boundary_2",
        type: "line",
        source: "borders",
        layout: { "line-cap": "round", "line-join": "round" },
        paint: {
          "line-color": "hsl(0,0%,70%)",
          "line-opacity": ["interpolate", ["linear"], ["zoom"], 0, 0.4, 4, 1],
          "line-width": ["interpolate", ["linear"], ["zoom"], 3, 1, 5, 1.8, 12, 3],
        },
      },
      labelLayer(3, 2, [3, 9, 7, 17]),
      labelLayer(2, 0, [2, 9, 5, 17]),
      labelLayer(1, 0, [1, 9, 4, 17]),
    ],
  };
}

function labelLayer(rank: number, minzoom: number, sizeStops: readonly number[]) {
  return {
    id: `label_country_${rank}`,
    type: "symbol",
    source: "labels",
    minzoom,
    maxzoom: 9,
    filter: ["==", ["get", "rank"], rank],
    layout: {
      "text-field": ["get", "name"],
      "text-font": ["Noto Sans Bold"],
      "text-max-width": 6.25,
      "text-size": ["interpolate", ["linear"], ["zoom"], ...sizeStops],
    },
    paint: {
      "text-color": LABEL_COLOR,
      "text-halo-blur": 1,
      "text-halo-color": "#fff",
      "text-halo-width": 1,
    },
  };
}
