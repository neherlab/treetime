# Toast checks of the upload tests pass without rendered portals

The upload tests check that no error toast appears with `queryByText(...)` equal to `null`. Toasts render through a Base UI portal, and in these tests the portal never renders, so the checks pass whether or not the code shows the toast. The "Example files" dialog, also in a portal, has no test that clicks through it.

## Evidence

- Checks for the absence of the toast "An earlier upload cannot be removed" [packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L154](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L154), [#L212](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L212)
- Toasts render in `ToastPortal` [packages/app-ui/src/ui/toast.tsx#L23-L25](../../packages/app-ui/src/ui/toast.tsx#L23-L25)
- Base UI chooses its layout effect when its module loads: `useIsoLayoutEffect` of `@base-ui/utils` is `React.useLayoutEffect` when `document` exists and a no-op otherwise. The portal node of `FloatingPortalLite` is resolved in that effect. The upload tests import the components statically and register happy-dom in `beforeAll` (`registerDom()`), after the imports have loaded Base UI, so every Base UI portal of these tests stays empty. A toast added under this setup does not reach the document
- `packages/app-ui/src/auspice/__tests__/test_notifications.test.tsx` registers happy-dom first and imports the toast components dynamically afterwards; there the toast text is found

## Fix direction

Load Base UI after the DOM in every test that renders Base UI components: register happy-dom before the components load (for example in a module that the test imports first, or with dynamic imports as in `test_notifications.test.tsx`), then pair every absence check with a test where the toast does appear, and add a click test of the "Example files" dialog.

## Validation

- Each toast absence check fails when the toast is forced to appear
