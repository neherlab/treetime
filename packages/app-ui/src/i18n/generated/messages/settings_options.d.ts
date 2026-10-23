export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_OptionsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_OptionsInputs */
/**
* | output |
* | --- |
* | "Values" |
*
* @param {Settings_OptionsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_options: ((inputs?: Settings_OptionsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_OptionsInputs, {
    locale?: "en" | "de";
}, {}>;
