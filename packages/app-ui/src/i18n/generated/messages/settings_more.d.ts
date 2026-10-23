export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_MoreInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_MoreInputs */
/**
* | output |
* | --- |
* | "More" |
*
* @param {Settings_MoreInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_more: ((inputs?: Settings_MoreInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_MoreInputs, {
    locale?: "en" | "de";
}, {}>;
