export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Settings_DescriptionInputs = {
    command: NonNullable<unknown>;
    title: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "Settings of treetime {command} ({title}), with their flags, config keys, and defaults." |
*
* @param {Docs_Settings_DescriptionInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_settings_description: ((inputs: Docs_Settings_DescriptionInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Settings_DescriptionInputs, {
    locale?: "en" | "de";
}, {}>;
