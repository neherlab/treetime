import { useState, useMemo } from "react";
import { Download } from "lucide-react";
import { clsx } from "clsx";
import { useAppStore } from "../store/app-store";
import type { CommandName } from "../types";

interface TabDef {
  key: string;
  label: string;
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

export function ResultsPanel() {
  const activeCommand = useAppStore((s) => s.activeCommand);
  const tabs = useMemo(() => COMMAND_TABS[activeCommand], [activeCommand]);
  const [activeTab, setActiveTab] = useState(tabs[0].key);

  const validTab = useMemo(() => (tabs.some((t) => t.key === activeTab) ? activeTab : tabs[0].key), [tabs, activeTab]);

  return (
    <div className="flex flex-1 flex-col overflow-hidden">
      <div className="flex items-center justify-between border-b border-gray-200 px-4 dark:border-gray-700">
        <div className="flex">
          {tabs.map((tab) => (
            <button
              key={tab.key}
              type="button"
              onClick={() => setActiveTab(tab.key)}
              className={clsx(
                "border-b-2 px-3 py-2 text-xs font-medium transition-colors",
                validTab === tab.key
                  ? "border-[var(--color-accent)] text-[var(--color-accent)]"
                  : "border-transparent text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300",
              )}
            >
              {tab.label}
            </button>
          ))}
        </div>
        <button
          type="button"
          className="flex items-center gap-1 rounded-md border border-gray-200 px-2 py-1 text-xs text-gray-500 hover:bg-gray-50 dark:border-gray-700 dark:text-gray-400 dark:hover:bg-gray-800"
        >
          <Download size={12} />
          Export
        </button>
      </div>

      <div className="flex-1 overflow-auto p-4">
        <TabContent command={activeCommand} tab={validTab} />
      </div>
    </div>
  );
}

function TabContent({ command, tab }: { command: CommandName; tab: string }) {
  if (tab === "tree") return <MockTree />;
  if (tab === "regression") return <MockRegressionPlot />;
  if (tab === "table" || tab === "traits" || tab === "confidence") return <MockDataTable command={command} tab={tab} />;
  if (tab === "sequences") return <MockFastaViewer />;
  if (tab === "model") return <MockModelPanel command={command} />;
  if (tab === "summary") return <MockPruneSummary />;
  if (tab === "auspice") return <MockAuspiceLink />;
  return <div className="text-sm text-gray-500">Unknown tab</div>;
}

function MockTree() {
  return (
    <div className="flex flex-col items-center justify-center rounded-lg border-2 border-dashed border-gray-200 bg-gray-50 p-8 dark:border-gray-700 dark:bg-gray-900">
      <div className="mb-4 font-mono text-xs text-gray-400 dark:text-gray-600">
        <pre className="leading-relaxed">{MOCK_TREE_ASCII}</pre>
      </div>
      <div className="flex gap-2">
        {["Rectangular", "Radial", "Clock"].map((layout) => (
          <button
            key={layout}
            type="button"
            className="rounded-md border border-gray-200 px-2.5 py-1 text-xs text-gray-500 hover:bg-white dark:border-gray-700 dark:text-gray-400 dark:hover:bg-gray-800"
          >
            {layout}
          </button>
        ))}
      </div>
      <p className="mt-3 text-xs text-gray-400">Tree visualization placeholder</p>
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
      <div className="rounded-lg border border-gray-200 bg-white p-4 dark:border-gray-700 dark:bg-gray-900">
        <svg viewBox="0 0 600 300" className="w-full">
          <line x1="60" y1="260" x2="570" y2="260" stroke="currentColor" className="text-gray-300 dark:text-gray-600" />
          <line x1="60" y1="20" x2="60" y2="260" stroke="currentColor" className="text-gray-300 dark:text-gray-600" />

          <text x="315" y="290" textAnchor="middle" className="fill-gray-500 text-[10px]">
            Sampling date
          </text>
          <text
            x="15"
            y="140"
            textAnchor="middle"
            transform="rotate(-90, 15, 140)"
            className="fill-gray-500 text-[10px]"
          >
            Root-to-tip divergence
          </text>

          <line
            x1="80"
            y1="240"
            x2="555"
            y2="40"
            stroke="currentColor"
            className="text-[var(--color-accent)]"
            strokeWidth="1.5"
            strokeDasharray="4 2"
          />

          {points.map((p) => (
            <circle
              key={p.name}
              cx={60 + ((p.x - 1999) / 15) * 510}
              cy={260 - (p.y / 0.05) * 240}
              r={p.outlier ? 5 : 3.5}
              className={p.outlier ? "fill-red-400" : "fill-[var(--color-accent)]"}
              opacity={0.7}
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
    <div className="rounded-lg border border-gray-200 bg-white px-3 py-2 dark:border-gray-700 dark:bg-gray-900">
      <div className="text-xs text-gray-500 dark:text-gray-400">{label}</div>
      <div className="text-sm font-semibold text-gray-900 dark:text-gray-100">{value}</div>
      {unit && <div className="text-xs text-gray-400">{unit}</div>}
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

function MockDataTable({ command, tab }: { command: CommandName; tab: string }) {
  const { columns, rows } = useMemo(() => getMockTableData(command, tab), [command, tab]);

  return (
    <div className="overflow-auto rounded-lg border border-gray-200 dark:border-gray-700">
      <table className="w-full text-left text-xs">
        <thead>
          <tr className="border-b border-gray-200 bg-gray-50 dark:border-gray-700 dark:bg-gray-800">
            {columns.map((col) => (
              <th key={col} className="px-3 py-2 font-medium text-gray-500 dark:text-gray-400">
                {col}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => (
            <tr
              key={row[0]}
              className="border-b border-gray-100 hover:bg-gray-50 dark:border-gray-800 dark:hover:bg-gray-800/50"
            >
              {row.map((cell, j) => (
                <td key={columns[j]} className="px-3 py-1.5 text-gray-700 dark:text-gray-300">
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

function getMockTableData(command: CommandName, tab: string): { columns: string[]; rows: string[][] } {
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
    <div className="overflow-auto rounded-lg border border-gray-200 font-mono text-xs dark:border-gray-700">
      {MOCK_FASTA.map((entry) => (
        <div key={entry.name} className="border-b border-gray-100 p-2 dark:border-gray-800">
          <div className="text-[var(--color-accent)]">&gt;{entry.name}</div>
          <div className="break-all text-gray-600 dark:text-gray-400">{entry.seq}</div>
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

function MockModelPanel({ command }: { command: CommandName }) {
  const showGtr = command !== "clock";

  return (
    <div className="space-y-4">
      {showGtr && (
        <div className="rounded-lg border border-gray-200 bg-white p-4 dark:border-gray-700 dark:bg-gray-900">
          <h4 className="mb-2 text-xs font-semibold text-gray-500 dark:text-gray-400">GTR Model: Inferred</h4>
          <div className="mb-3">
            <div className="mb-1 text-xs text-gray-500">Rate matrix</div>
            <table className="font-mono text-xs">
              <thead>
                <tr>
                  <th className="w-8" />
                  {NUC_LABELS.map((n) => (
                    <th key={n} className="w-16 px-2 text-center text-gray-500">
                      {n}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {NUC_LABELS.map((row, i) => (
                  <tr key={row}>
                    <td className="pr-2 text-gray-500">{row}</td>
                    {GTR_RATE_MATRIX[i].map((val, colIdx) => (
                      <td key={NUC_LABELS[colIdx]} className="px-2 text-center text-gray-700 dark:text-gray-300">
                        {i === colIdx ? "-" : val.toFixed(2)}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div>
            <div className="mb-1 text-xs text-gray-500">Equilibrium frequencies</div>
            <div className="flex gap-4 font-mono text-xs text-gray-700 dark:text-gray-300">
              {GTR_FREQUENCIES.map(({ nuc, freq }) => (
                <div key={nuc} className="flex items-center gap-1.5">
                  <span className="text-gray-500">{nuc}:</span>
                  <div className="h-3 w-16 rounded-sm bg-gray-200 dark:bg-gray-700">
                    <div className="h-full rounded-sm bg-[var(--color-accent)]" style={{ width: `${freq * 100}%` }} />
                  </div>
                  <span>{freq.toFixed(2)}</span>
                </div>
              ))}
            </div>
          </div>
        </div>
      )}

      {(command === "clock" || command === "timetree") && (
        <div className="rounded-lg border border-gray-200 bg-white p-4 dark:border-gray-700 dark:bg-gray-900">
          <h4 className="mb-2 text-xs font-semibold text-gray-500 dark:text-gray-400">Clock Model</h4>
          <div className="space-y-1 text-xs">
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

function KV({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex items-center gap-2">
      <span className="w-24 text-gray-500 dark:text-gray-400">{label}</span>
      <span className="font-mono text-gray-700 dark:text-gray-300">{value}</span>
    </div>
  );
}

function MockPruneSummary() {
  return (
    <div className="rounded-lg border border-gray-200 bg-white p-4 dark:border-gray-700 dark:bg-gray-900">
      <h4 className="mb-3 text-xs font-semibold text-gray-500 dark:text-gray-400">Pruning Summary</h4>
      <div className="space-y-1 text-xs">
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
    <div className="flex flex-col items-center gap-3 rounded-lg border-2 border-dashed border-gray-200 p-8 dark:border-gray-700">
      <p className="text-sm text-gray-600 dark:text-gray-400">Open the Auspice v2 JSON in Nextstrain's tree viewer</p>
      <button
        type="button"
        className="rounded-lg bg-[var(--color-accent)] px-4 py-2 text-sm font-medium text-white hover:bg-[var(--color-accent-hover)]"
      >
        Open in Auspice
      </button>
      <p className="text-xs text-gray-400">auspice_tree.json will be generated after a real run</p>
    </div>
  );
}
