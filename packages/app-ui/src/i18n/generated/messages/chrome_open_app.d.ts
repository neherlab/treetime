export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Open_AppInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_Open_AppInputs */
/**
* | output |
* | --- |
* | "Open TreeTime" |
*
* @param {Chrome_Open_AppInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_open_app: ((inputs?: Chrome_Open_AppInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Open_AppInputs, {
    locale?: "en" | "de";
}, {}>;
