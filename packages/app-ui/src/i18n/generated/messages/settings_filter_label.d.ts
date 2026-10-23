export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Filter_LabelInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Filter_LabelInputs */
/**
* | output |
* | --- |
* | "Filter settings" |
*
* @param {Settings_Filter_LabelInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_filter_label: ((inputs?: Settings_Filter_LabelInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_LabelInputs, {
    locale?: "en" | "de";
}, {}>;
