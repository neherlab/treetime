export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_InputsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_InputsInputs */
/**
* | output |
* | --- |
* | "Inputs" |
*
* @param {Settings_InputsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_inputs: ((inputs?: Settings_InputsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_InputsInputs, {
    locale?: "en" | "de";
}, {}>;
