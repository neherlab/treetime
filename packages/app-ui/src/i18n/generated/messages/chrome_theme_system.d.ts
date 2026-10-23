export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Theme_SystemInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_Theme_SystemInputs */
/**
* | output |
* | --- |
* | "Theme follows the system" |
*
* @param {Chrome_Theme_SystemInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_theme_system: ((inputs?: Chrome_Theme_SystemInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_SystemInputs, {
    locale?: "en" | "de";
}, {}>;
