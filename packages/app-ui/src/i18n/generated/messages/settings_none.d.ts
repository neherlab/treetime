export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_NoneInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_NoneInputs */
/**
* | output |
* | --- |
* | "none" |
*
* @param {Settings_NoneInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_none: ((inputs?: Settings_NoneInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_NoneInputs, {
    locale?: "en" | "de";
}, {}>;
