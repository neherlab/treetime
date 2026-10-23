export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Theme_DarkInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_Theme_DarkInputs */
/**
* | output |
* | --- |
* | "Dark theme" |
*
* @param {Chrome_Theme_DarkInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_theme_dark: ((inputs?: Chrome_Theme_DarkInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_DarkInputs, {
    locale?: "en" | "de";
}, {}>;
