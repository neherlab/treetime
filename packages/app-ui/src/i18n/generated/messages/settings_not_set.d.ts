export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Not_SetInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Not_SetInputs */
/**
* | output |
* | --- |
* | "not set" |
*
* @param {Settings_Not_SetInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_not_set: ((inputs?: Settings_Not_SetInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Not_SetInputs, {
    locale?: "en" | "de";
}, {}>;
