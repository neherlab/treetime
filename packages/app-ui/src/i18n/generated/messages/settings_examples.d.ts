export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_ExamplesInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_ExamplesInputs */
/**
* | output |
* | --- |
* | "Examples" |
*
* @param {Settings_ExamplesInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_examples: ((inputs?: Settings_ExamplesInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_ExamplesInputs, {
    locale?: "en" | "de";
}, {}>;
