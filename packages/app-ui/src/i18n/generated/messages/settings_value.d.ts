export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_ValueInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_ValueInputs */
/**
* | output |
* | --- |
* | "Value" |
*
* @param {Settings_ValueInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_value: ((inputs?: Settings_ValueInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_ValueInputs, {
    locale?: "en" | "de";
}, {}>;
