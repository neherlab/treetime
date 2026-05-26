import { Header } from "./components/Header";
import { CommandNav } from "./components/CommandNav";
import { Workspace } from "./components/Workspace";
import { useElectronThemeSync } from "./hooks/useElectronThemeSync";

export function App() {
  useElectronThemeSync();

  return (
    <div className="flex h-screen flex-col overflow-hidden bg-white dark:bg-gray-950">
      <Header />
      <div className="flex flex-1 overflow-hidden">
        <CommandNav />
        <Workspace />
      </div>
    </div>
  );
}
