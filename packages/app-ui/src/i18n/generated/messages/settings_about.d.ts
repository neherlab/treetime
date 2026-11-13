export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_AboutInputs = {
    label: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "About {label}" |
*
* @param {Settings_AboutInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_about: ((inputs: Settings_AboutInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_AboutInputs, {
    locale?: "en" | "de";
}, {}>;
