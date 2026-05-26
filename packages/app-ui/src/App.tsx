import { useVersion } from "./hooks";

export function App() {
  const { data: version, isLoading, error } = useVersion();

  return (
    <div style={{ fontFamily: "system-ui, sans-serif", padding: "2rem" }}>
      <h1>TreeTime</h1>
      <p>Phylodynamic inference</p>
      {isLoading && <p>Loading...</p>}
      {error && <p style={{ color: "red" }}>Failed to connect to server: {error.message}</p>}
      {version && <p>v{version.version}</p>}
    </div>
  );
}
