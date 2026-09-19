import { RouterProvider } from "@tanstack/react-router";

import "./ui/fonts";
import { router } from "./router";

export function App() {
  return <RouterProvider router={router} />;
}
