export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Theme_LightInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_Theme_LightInputs */
/**
* | output |
* | --- |
* | "Light theme" |
*
* @param {Chrome_Theme_LightInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_theme_light: ((inputs?: Chrome_Theme_LightInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_LightInputs, {
    locale?: "en" | "de";
}, {}>;
