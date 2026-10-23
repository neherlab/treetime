export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Filter_PlaceholderInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Filter_PlaceholderInputs */
/**
* | output |
* | --- |
* | "Label, flag, config key, or help text" |
*
* @param {Settings_Filter_PlaceholderInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_filter_placeholder: ((inputs?: Settings_Filter_PlaceholderInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_PlaceholderInputs, {
    locale?: "en" | "de";
}, {}>;
