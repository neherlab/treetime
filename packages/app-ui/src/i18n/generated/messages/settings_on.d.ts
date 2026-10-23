export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_OnInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_OnInputs */
/**
* | output |
* | --- |
* | "on" |
*
* @param {Settings_OnInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_on: ((inputs?: Settings_OnInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_OnInputs, {
    locale?: "en" | "de";
}, {}>;
