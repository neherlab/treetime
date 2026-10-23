export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Need_RequiredInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Need_RequiredInputs */
/**
* | output |
* | --- |
* | "required" |
*
* @param {Settings_Need_RequiredInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_need_required: ((inputs?: Settings_Need_RequiredInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Need_RequiredInputs, {
    locale?: "en" | "de";
}, {}>;
