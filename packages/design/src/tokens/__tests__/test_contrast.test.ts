import { contrastRatio } from "../color";
import { CONTRAST_REQUIREMENTS, resolveColor } from "../tokens";

const THEMES = ["light", "dark"] as const;

describe("token contrast", () => {
  test.each(
    CONTRAST_REQUIREMENTS.flatMap(({ foreground, background, minimum }) =>
      THEMES.map((theme) => ({ foreground, background, minimum, theme })),
    ),
  )(
    "$foreground on $background in the $theme theme reaches $minimum:1",
    ({ foreground, background, minimum, theme }) => {
      const ratio = contrastRatio(resolveColor(foreground)[theme], resolveColor(background)[theme]);

      expect(ratio).toBeGreaterThanOrEqual(minimum);
    },
  );
});
