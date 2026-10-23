export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_FieldsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_FieldsInputs */
/**
* | output |
* | --- |
* | "Fields" |
*
* @param {Settings_FieldsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_fields: ((inputs?: Settings_FieldsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_FieldsInputs, {
    locale?: "en" | "de";
}, {}>;
