export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Need_OptionalInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Need_OptionalInputs */
/**
* | output |
* | --- |
* | "optional" |
*
* @param {Settings_Need_OptionalInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_need_optional: ((inputs?: Settings_Need_OptionalInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Need_OptionalInputs, {
    locale?: "en" | "de";
}, {}>;
