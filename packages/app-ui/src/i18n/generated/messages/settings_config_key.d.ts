export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Config_KeyInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Config_KeyInputs */
/**
* | output |
* | --- |
* | "Config key" |
*
* @param {Settings_Config_KeyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_config_key: ((inputs?: Settings_Config_KeyInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Config_KeyInputs, {
    locale?: "en" | "de";
}, {}>;
