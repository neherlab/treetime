# Toast checks of the upload tests pass without rendered portals

The upload tests check that no error toast appears with `queryByText(...)` equal to `null`. Toasts render through a Base UI portal, and in these tests the portal never renders, so the checks pass whether or not the code shows the toast. The "Example files" dialog, also in a portal, has no test that clicks through it.

## Evidence

- Checks for the absence of the toast "An earlier upload cannot be removed" [`packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L150`](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L150), [#L208](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L208)
- Toasts render in `ToastPortal` [packages/design/src/ui/toast.tsx#L23-L25](../../packages/design/src/ui/toast.tsx#L23-L25)
- Base UI chooses its layout effect when its module loads: `useIsoLayoutEffect` of `@base-ui/utils` is `React.useLayoutEffect` when `document` exists and a no-op otherwise. The portal node of `FloatingPortalLite` is resolved in that effect. The upload tests import the components statically and register happy-dom in `beforeAll` (`registerDom()`), after the imports have loaded Base UI, so every Base UI portal of these tests stays empty. A toast added under this setup does not reach the document
- `packages/app-ui/src/auspice/__tests__/test_notifications.test.tsx` registers happy-dom first and imports the toast components dynamically afterwards; there the toast text is found

## Fix direction

Load Base UI after the DOM in every test that renders Base UI components: register happy-dom before the components load (for example in a module that the test imports first, or with dynamic imports as in `test_notifications.test.tsx`), then pair every absence check with a test where the toast does appear, and add a click test of the "Example files" dialog.

## Validation

- Each toast absence check fails when the toast is forced to appear

## Proposed implementation plan

1. Make happy-dom exist before Base UI loads in the upload tests: move `registerDom()` of `packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx` (now in `beforeAll`, line 39) into a setup module that the test imports first, or load `Toaster` and `InputsSection` with dynamic imports after `registerDom()`, as `packages/app-ui/src/auspice/__tests__/test_notifications.test.tsx` (lines 11, 70-71) does. Prefer one shared setup module in `packages/app-ui/src/__tests__/` (beside `dom.ts`) so that every Base UI test gets the same order.
2. Add a positive sibling case for each absence check at lines 150 and 208: a server fake whose release fails, and an assertion that the toast "An earlier upload cannot be removed" appears. Oracle: the error path of the upload code that shows this toast.
3. Add a test that opens the "Example files" dialog with `userEvent`, picks an example, and asserts the resulting upload entry by role or label.
4. Validate: `just test-ts`; see red by forcing the toast to always appear and confirming that the checks at lines 150 and 208 fail, then revert.
