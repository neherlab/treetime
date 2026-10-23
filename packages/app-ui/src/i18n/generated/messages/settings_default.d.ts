export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_DefaultInputs = {
    value: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "Default: {value}" |
*
* @param {Settings_DefaultInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_default: ((inputs: Settings_DefaultInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_DefaultInputs, {
    locale?: "en" | "de";
}, {}>;
