import { Download } from "lucide-react";
import { useTheme } from "next-themes";
import { useCallback, useMemo, useState } from "react";

import { useActiveCommand } from "../hooks/useActiveCommand";
import type { CommandName } from "../types";
import { Button, Menu, cn } from "../ui";

interface TabDef {
  key: string;
  label: string;
}

interface TableData {
  columns: string[];
  rows: string[][];
}

const COMMAND_TABS: Record<CommandName, TabDef[]> = {
  timetree: [
    { key: "tree", label: "Tree" },
    { key: "model", label: "Model" },
    { key: "confidence", label: "Confidence" },
    { key: "auspice", label: "Auspice" },
  ],
  ancestral: [
    { key: "tree", label: "Tree" },
    { key: "sequences", label: "Sequences" },
    { key: "model", label: "Model" },
  ],
  clock: [
    { key: "regression", label: "Regression" },
    { key: "tree", label: "Tree" },
    { key: "table", label: "Table" },
    { key: "model", label: "Model" },
  ],
  mugration: [
    { key: "tree", label: "Tree" },
    { key: "traits", label: "Traits" },
    { key: "model", label: "Model" },
  ],
  optimize: [
    { key: "tree", label: "Tree" },
    { key: "model", label: "Model" },
  ],
  prune: [
    { key: "tree", label: "Tree" },
    { key: "summary", label: "Summary" },
  ],
};

const EXPORT_TARGETS = ["Auspice JSON", "Newick tree", "Node data (CSV)"];

export function ResultsPanel() {
  const activeCommand = useActiveCommand();
  const tabs = useMemo(() => COMMAND_TABS[activeCommand], [activeCommand]);
  const fallbackTab = tabs[0]?.key ?? "tree";
  const [activeTab, setActiveTab] = useState(fallbackTab);

  const validTab = useMemo(
    () => (tabs.some((t) => t.key === activeTab) ? activeTab : fallbackTab),
    [tabs, activeTab, fallbackTab],
  );

  return (
    <div className="flex flex-1 flex-col overflow-hidden">
      <div className="border-line flex items-center justify-between border-b px-4">
        <div className="flex">
          {tabs.map((tab) => (
            <TabButton key={tab.key} tab={tab} active={validTab === tab.key} onSelect={setActiveTab} />
          ))}
        </div>
        <Menu.Root>
          <Menu.Trigger
            render={
              <Button variant="outline" size="sm">
                <Download size={12} />
                Export
              </Button>
            }
          />
          <Menu.Popup>
            {EXPORT_TARGETS.map((target) => (
              <Menu.Item key={target}>{target}</Menu.Item>
            ))}
          </Menu.Popup>
        </Menu.Root>
      </div>

      <div className="flex-1 overflow-auto p-4">
        <TabContent command={activeCommand} tab={validTab} />
      </div>
    </div>
  );
}

function TabButton({ tab, active, onSelect }: { tab: TabDef; active: boolean; onSelect: (key: string) => void }) {
  const handleClick = useCallback(() => onSelect(tab.key), [onSelect, tab.key]);

  return (
    <button
      type="button"
      onClick={handleClick}
      className={cn(
        "text-2xs border-b-2 px-3 py-2 font-medium transition-colors",
        active ? "border-accent text-accent" : "text-ink-muted hover:text-ink border-transparent",
      )}
    >
      {tab.label}
    </button>
  );
}

function TabContent({ command, tab }: { command: CommandName; tab: string }) {
  if (tab === "tree") return <MockTree />;

  if (tab === "regression") return <MockRegressionPlot />;

  if (tab === "table" || tab === "traits" || tab === "confidence") return <MockDataTable tab={tab} />;

  if (tab === "sequences") return <MockFastaViewer />;

  if (tab === "model") return <MockModelPanel command={command} />;

  if (tab === "summary") return <MockPruneSummary />;

  if (tab === "auspice") return <MockAuspiceLink />;

  return <div className="text-ink-muted text-sm">Unknown tab</div>;
}

function MockTree() {
  return (
    <div className="border-line bg-surface-1 flex flex-col items-center justify-center rounded-md border border-dashed p-8">
      <div className="text-2xs text-ink-faint mb-4 font-mono">
        <pre className="leading-relaxed">{MOCK_TREE_ASCII}</pre>
      </div>
      <div className="flex gap-2">
        {["Rectangular", "Radial", "Clock"].map((layout) => (
          <Button key={layout} variant="subtle" size="sm">
            {layout}
          </Button>
        ))}
      </div>
      <p className="text-2xs text-ink-faint mt-3">Tree visualization placeholder</p>
    </div>
  );
}

const MOCK_TREE_ASCII = `\
         +-- A/Indiana/03/2012
    +----+
    |    |    +-- A/Peru/PER247/2011
    |    +----+
    |         |    +-- A/Minab/797/2011
    |         +----+
    |              |    +-- A/Oregon/15/2009
    |              +----+
    |              |    +-- A/Hong_Kong/H090/2009
    |              |
    |              +-- A/Boston/57/2008
    |                  +-- A/DaNang/DN434/2008
    |                      +-- A/Managua/25/2007
    |                          +-- A/Mexico/InDRE940/2003
    |                              +-- A/New_York/182/2000
    |                              |   +-- A/Scotland/76/2003
    |                              |   +-- A/Denmark/107/2003
    |                              +-- A/Canterbury/58/2000
----+
    |    +-- A/Nebraska/15/2011
    +----+
    |    +-- A/Maryland/21/2011
    |
    +-- A/Maryland/03/2013
    +-- A/New_Hampshire/12/2012
    |
    +-- A/Hawaii/02/2013
    +-- A/Boston/DOA2_107/2012`;

function MockRegressionPlot() {
  const points = useMemo(() => MOCK_CLOCK_DATA, []);

  return (
    <div className="space-y-3">
      <div className="border-line bg-surface-1 rounded-md border p-4">
        <svg viewBox="0 0 600 300" className="w-full">
          <line x1="60" y1="260" x2="570" y2="260" stroke="currentColor" className="text-line-strong" />
          <line x1="60" y1="20" x2="60" y2="260" stroke="currentColor" className="text-line-strong" />

          <text x="315" y="290" textAnchor="middle" className="fill-ink-muted text-[10px]">
            Sampling date
          </text>
          <text
            x="15"
            y="140"
            textAnchor="middle"
            transform="rotate(-90, 15, 140)"
            className="fill-ink-muted text-[10px]"
          >
            Root-to-tip divergence
          </text>

          <line
            x1="80"
            y1="240"
            x2="555"
            y2="40"
            stroke="currentColor"
            className="text-accent"
            strokeWidth="1.5"
            strokeDasharray="4 2"
          />

          {points.map((p) => (
            <circle
              key={p.name}
              cx={60 + ((p.x - 1999) / 15) * 510}
              cy={260 - (p.y / 0.05) * 240}
              r={p.outlier ? 5 : 3.5}
              className={p.outlier ? "fill-signal-outlier" : "fill-accent"}
              opacity={0.75}
            />
          ))}
        </svg>
      </div>

      <div className="grid grid-cols-3 gap-3">
        <StatCard label="Clock rate" value="3.3e-3" unit="subs/site/year" />
        <StatCard label="R-squared" value="0.97" />
        <StatCard label="Intercept" value="-6.57" />
      </div>
    </div>
  );
}

function StatCard({ label, value, unit }: { label: string; value: string; unit?: string }) {
  return (
    <div className="border-line bg-surface-1 rounded-md border px-3 py-2">
      <div className="text-2xs text-ink-muted">{label}</div>
      <div className="text-ink font-mono text-sm font-medium">{value}</div>
      {unit && <div className="text-2xs text-ink-faint">{unit}</div>}
    </div>
  );
}

const MOCK_CLOCK_DATA = [
  { x: 2000.134, y: 0.0021, name: "A/New_York/182/2000", outlier: false },
  { x: 2000.682, y: 0.0072, name: "A/Canterbury/58/2000", outlier: false },
  { x: 2003.003, y: 0.0098, name: "A/Mexico/InDRE940/2003", outlier: false },
  { x: 2003.003, y: 0.0253, name: "A/Denmark/107/2003", outlier: true },
  { x: 2003.841, y: 0.0251, name: "A/Scotland/76/2003", outlier: true },
  { x: 2007.487, y: 0.0193, name: "A/Managua/25/2007", outlier: false },
  { x: 2008.151, y: 0.0206, name: "A/Boston/57/2008", outlier: false },
  { x: 2008.865, y: 0.0249, name: "A/DaNang/DN434/2008", outlier: false },
  { x: 2009.482, y: 0.027, name: "A/Oregon/15/2009", outlier: false },
  { x: 2009.523, y: 0.0284, name: "A/Hong_Kong/H090/2009", outlier: false },
  { x: 2011.652, y: 0.0338, name: "A/Peru/PER247/2011", outlier: false },
  { x: 2011.956, y: 0.0371, name: "A/Nebraska/15/2011", outlier: false },
  { x: 2011.98, y: 0.0362, name: "A/Minab/797/2011", outlier: false },
  { x: 2011.986, y: 0.0378, name: "A/Maryland/21/2011", outlier: false },
  { x: 2012.257, y: 0.0385, name: "A/Indiana/03/2012", outlier: false },
  { x: 2012.838, y: 0.0399, name: "A/Boston/DOA2_107/2012", outlier: false },
  { x: 2012.857, y: 0.0427, name: "A/New_Hampshire/12/2012", outlier: false },
  { x: 2013.112, y: 0.0456, name: "A/Maryland/03/2013", outlier: false },
  { x: 2013.405, y: 0.0441, name: "A/Hawaii/02/2013", outlier: false },
];

function MockDataTable({ tab }: { tab: string }) {
  const { columns, rows } = useMemo(() => getMockTableData(tab), [tab]);

  return (
    <div className="border-line overflow-auto rounded-md border">
      <table className="text-2xs w-full text-left">
        <thead>
          <tr className="border-line bg-surface-2 border-b">
            {columns.map((col) => (
              <th key={col} className="text-ink-muted px-3 py-2 font-medium">
                {col}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => (
            <tr key={row[0]} className="border-line/60 hover:bg-surface-1 border-b">
              {row.map((cell, j) => (
                <td key={columns[j] ?? j} className="text-ink px-3 py-1.5 font-mono">
                  {cell}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function getMockTableData(tab: string): TableData {
  if (tab === "table") {
    return {
      columns: ["Name", "Divergence", "Date", "Predicted", "Deviation", "Outlier"],
      rows: MOCK_CLOCK_DATA.map((p) => [
        p.name,
        p.y.toFixed(4),
        p.x.toFixed(1),
        (p.x + 0.1).toFixed(1),
        ((p.y - 0.005 * (p.x - 2018)) * 100).toFixed(2),
        p.outlier ? "Yes" : "",
      ]),
    };
  }

  if (tab === "traits") {
    return {
      columns: ["Node", "Trait", "Confidence"],
      rows: [
        ["A/Hawaii/02/2013", "USA", "0.98"],
        ["A/Boston/DOA2_107/2012", "USA", "0.97"],
        ["A/Oregon/15/2009", "USA", "0.96"],
        ["A/Hong_Kong/H090/2009", "Hong Kong", "0.94"],
        ["A/Canterbury/58/2000", "New Zealand", "0.92"],
        ["A/Managua/25/2007", "Nicaragua", "0.91"],
        ["A/DaNang/DN434/2008", "Viet Nam", "0.89"],
        ["A/Peru/PER247/2011", "Peru", "0.93"],
        ["A/Minab/797/2011", "Iran", "0.88"],
        ["A/Denmark/107/2003", "Denmark", "0.95"],
        ["A/Scotland/76/2003", "United Kingdom", "0.90"],
        ["NODE_0001", "USA", "0.72"],
        ["NODE_0002", "USA", "0.65"],
        ["NODE_0003", "Asia", "0.58"],
      ],
    };
  }

  if (tab === "confidence") {
    return {
      columns: ["Node", "Date", "Lower CI", "Upper CI"],
      rows: [
        ["NODE_0001", "2011.4", "2010.9", "2011.8"],
        ["NODE_0002", "2009.1", "2008.3", "2009.7"],
        ["NODE_0003", "2006.8", "2005.9", "2007.5"],
        ["NODE_0004", "2003.2", "2002.1", "2004.0"],
        ["NODE_0005", "2000.5", "1999.2", "2001.4"],
        ["root", "1998.7", "1997.1", "1999.8"],
      ],
    };
  }

  return { columns: ["Column"], rows: [["No data"]] };
}

function MockFastaViewer() {
  return (
    <div className="border-line text-2xs overflow-auto rounded-md border font-mono">
      {MOCK_FASTA.map((entry) => (
        <div key={entry.name} className="border-line/60 border-b p-2">
          <div className="text-accent">&gt;{entry.name}</div>
          <div className="text-ink-muted break-all">{entry.seq}</div>
        </div>
      ))}
    </div>
  );
}

const MOCK_FASTA = [
  {
    name: "A/Hawaii/02/2013",
    seq: "ATGAATCCAAATCAAAAGATAATAACAATTGGCTCTGTTTCTCTCACCATTTCCACAGTATGCTTCTTCATGCAAATTGC...",
  },
  {
    name: "A/Indiana/03/2012",
    seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGC...",
  },
  {
    name: "A/Oregon/15/2009",
    seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGC...",
  },
  {
    name: "A/New_York/182/2000",
    seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTGCCACAATATGCTTCCTTATGCAAATTGC...",
  },
  { name: "NODE_0001", seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGC..." },
  { name: "NODE_0002", seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTBCCACAATATGCTTCYTCATGCAAATTGC..." },
  { name: "root", seq: "ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTBCCACAATATGCTTCYTCWTGCAAATTGC..." },
];

const NUC_LABELS = ["A", "C", "G", "T"];

const GTR_RATE_MATRIX = [
  [0, 0.94, 2.41, 0.52],
  [0.94, 0, 0.48, 2.68],
  [2.41, 0.48, 0, 0.87],
  [0.52, 2.68, 0.87, 0],
];

const GTR_FREQUENCIES = [
  { nuc: "A", freq: 0.334 },
  { nuc: "C", freq: 0.198 },
  { nuc: "G", freq: 0.223 },
  { nuc: "T", freq: 0.245 },
];

const FREQ_HUES = [200, 150, 70, 25];

function MockModelPanel({ command }: { command: CommandName }) {
  const showGtr = command !== "clock";
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === "dark";

  const maxRate = useMemo(() => Math.max(...GTR_RATE_MATRIX.flatMap((row, i) => row.filter((_, j) => i !== j))), []);

  return (
    <div className="space-y-4">
      {showGtr && (
        <div className="border-line bg-surface-1 rounded-md border p-4">
          <h4 className="text-2xs text-ink-muted mb-3 font-semibold">GTR model (inferred)</h4>
          <div className="mb-4">
            <div className="text-2xs text-ink-muted mb-1.5">Rate matrix</div>
            <RateMatrix maxRate={maxRate} isDark={isDark} />
          </div>
          <div>
            <div className="text-2xs text-ink-muted mb-1.5">Equilibrium frequencies</div>
            <FrequencyBar isDark={isDark} />
          </div>
        </div>
      )}

      {(command === "clock" || command === "timetree") && (
        <div className="border-line bg-surface-1 rounded-md border p-4">
          <h4 className="text-2xs text-ink-muted mb-2 font-semibold">Clock model</h4>
          <div className="text-2xs space-y-1">
            <KV label="Clock rate" value="3.3e-3 subs/site/year" />
            <KV label="Intercept" value="-6.57" />
            <KV label="R-squared" value="0.970" />
            <KV label="Chi-squared" value="8.7" />
          </div>
        </div>
      )}
    </div>
  );
}

function RateMatrix({ maxRate, isDark }: { maxRate: number; isDark: boolean }) {
  const cell = 34;
  const pad = 18;
  const size = pad + cell * NUC_LABELS.length;

  return (
    <svg viewBox={`0 0 ${size} ${size}`} className="w-full max-w-[220px] font-mono">
      {NUC_LABELS.map((label, j) => (
        <text
          key={`col-${label}`}
          x={pad + cell * j + cell / 2}
          y={pad - 6}
          textAnchor="middle"
          className="fill-ink-muted text-[9px]"
        >
          {label}
        </text>
      ))}
      {NUC_LABELS.map((label, i) => (
        <text
          key={`row-${label}`}
          x={pad - 6}
          y={pad + cell * i + cell / 2 + 3}
          textAnchor="end"
          className="fill-ink-muted text-[9px]"
        >
          {label}
        </text>
      ))}
      {NUC_LABELS.map((_, i) =>
        (GTR_RATE_MATRIX[i] ?? []).map((val, j) => {
          const isDiag = i === j;

          return (
            <g key={`${NUC_LABELS[i]}-${NUC_LABELS[j]}`}>
              <rect
                x={pad + cell * j + 1}
                y={pad + cell * i + 1}
                width={cell - 2}
                height={cell - 2}
                rx={3}
                fill={isDiag ? "var(--color-surface-2)" : rateHeatmapColor(val, maxRate, isDark)}
              />
              <text
                x={pad + cell * j + cell / 2}
                y={pad + cell * i + cell / 2 + 3}
                textAnchor="middle"
                className={isDiag ? "fill-ink-faint text-[9px]" : "fill-ink text-[9px]"}
              >
                {isDiag ? "-" : val.toFixed(2)}
              </text>
            </g>
          );
        }),
      )}
    </svg>
  );
}

function rateHeatmapColor(value: number, maxValue: number, isDark: boolean): string {
  const t = Math.min(value / maxValue, 1);
  const lightness = isDark ? 0.3 + 0.15 * (1 - t) : 0.92 - 0.15 * t;
  const chroma = 0.04 + 0.08 * t;
  const hue = 200 - 175 * t;

  return `oklch(${lightness.toFixed(3)} ${chroma.toFixed(3)} ${hue.toFixed(0)})`;
}

function FrequencyBar({ isDark }: { isDark: boolean }) {
  return (
    <svg viewBox="0 0 100 12" className="h-5 w-full" preserveAspectRatio="none">
      {GTR_FREQUENCIES.map(({ nuc, freq }, index) => {
        const x = GTR_FREQUENCIES.slice(0, index).reduce((sum, g) => sum + g.freq * 100, 0);
        const width = freq * 100;

        return (
          <g key={nuc}>
            <rect x={x} y={0} width={width} height={12} fill={freqColor(index, isDark)} />
            <text x={x + width / 2} y={8.5} textAnchor="middle" className="fill-ink text-[6px] font-medium">
              {nuc} {(freq * 100).toFixed(0)}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

function freqColor(index: number, isDark: boolean): string {
  const hue = FREQ_HUES[index] ?? 200;

  return isDark ? `oklch(0.5 0.08 ${hue})` : `oklch(0.82 0.08 ${hue})`;
}

function KV({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex items-center gap-2">
      <span className="text-ink-muted w-24">{label}</span>
      <span className="text-ink font-mono">{value}</span>
    </div>
  );
}

function MockPruneSummary() {
  return (
    <div className="border-line bg-surface-1 rounded-md border p-4">
      <h4 className="text-2xs text-ink-muted mb-3 font-semibold">Pruning summary</h4>
      <div className="text-2xs space-y-1">
        <KV label="Nodes removed" value="4" />
        <KV label="Branches removed" value="4" />
        <KV label="Nodes before" value="37" />
        <KV label="Nodes after" value="33" />
        <KV label="Leaves before" value="19" />
        <KV label="Leaves after" value="17" />
      </div>
    </div>
  );
}

function MockAuspiceLink() {
  return (
    <div className="border-line flex flex-col items-center gap-3 rounded-md border border-dashed p-8">
      <p className="text-ink-muted text-sm">Open the Auspice v2 JSON in Nextstrain's tree viewer</p>
      <Button variant="solid" size="md">
        Open in Auspice
      </Button>
      <p className="text-2xs text-ink-faint">auspice_tree.json is written after a real run</p>
    </div>
  );
}
