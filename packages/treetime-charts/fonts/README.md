# Chart fonts

`treetime-charts` embeds these fonts, draws chart text with them, and subsets them into PDF charts.

- **IBM Plex Sans** is the chart font: the interface typeface of the neherlab design language (`packages/design`), which the web app loads through `@fontsource-variable/ibm-plex-sans`. The charts need its static Regular and Bold instances, because the chart library and the PDF writer draw a variable font only at its default weight
- **Lato** draws the characters that IBM Plex Sans lacks, such as the superscript minus `⁻` (U+207B) of exponents like `10⁻³`. The SVG renderer and the PDF writer take a missing glyph from the next font of the font database, and `text_width()` measures each character in the font that draws it

## IBM Plex Sans

- Source: the release [`@ibm/plex-sans@1.1.0`](https://github.com/IBM/plex/releases/tag/%40ibm%2Fplex-sans%401.1.0) of [`IBM/plex`](https://github.com/IBM/plex), asset `ibm-plex-sans.zip` (sha256 `fb365d910566e6d199cc2c15579a7dd9a267128e18431a394ed81f1970c69200`), files `fonts/complete/ttf/IBMPlexSans-Regular.ttf`, `fonts/complete/ttf/IBMPlexSans-Bold.ttf`, and `LICENSE.txt`
- License: SIL Open Font License 1.1 with the Reserved Font Name "Plex" (`IBMPlexSans-OFL.txt`). The files are unmodified

## Lato

- Source: [`google/fonts`](https://github.com/google/fonts), folder [`ofl/lato/`](https://github.com/google/fonts/tree/5d3b76120a319730fda218cc7410174a462b32cb/ofl/lato) at commit `5d3b76120a319730fda218cc7410174a462b32cb`
- License: SIL Open Font License 1.1 with the Reserved Font Name "Lato" (`Lato-OFL.txt`). The files are unmodified

## Checksums

| File                      | sha256                                                             |
| ------------------------- | ------------------------------------------------------------------ |
| `IBMPlexSans-Regular.ttf` | `975dcda37d80f038dcd143c22e33ca2d97a0cc5a929aace1c749153b0fe1afa5` |
| `IBMPlexSans-Bold.ttf`    | `9e6c74a889a700d707613d24548fe4ffa6bc59559a0689d2cf9e133bdcdafb2f` |
| `IBMPlexSans-OFL.txt`     | `7e6b2818edbd8f6a01ae80641cc8f16a51080d08fb4e532be3a0b6f74adb07da` |
| `Lato-Regular.ttf`        | `d636e4683231f931eda222d588e944d082bfd3bdba02f928bee461c0f185b251` |
| `Lato-Bold.ttf`           | `8a0aace75d33794eece4b28187bfc1df0bbd2888b5d8a56e01788c8d65d16be1` |
| `Lato-OFL.txt`            | `74ba064d03f1f1c4a952da936c3eb71866c34404916734de3cae73b34357e59e` |
