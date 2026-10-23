export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_OffInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_OffInputs */
/**
* | output |
* | --- |
* | "off" |
*
* @param {Settings_OffInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_off: ((inputs?: Settings_OffInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_OffInputs, {
    locale?: "en" | "de";
}, {}>;
