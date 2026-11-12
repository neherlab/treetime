export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Command_LineInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_Command_LineInputs */
/**
* | output |
* | --- |
* | "Command line" |
*
* @param {Settings_Command_LineInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_command_line: ((inputs?: Settings_Command_LineInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Command_LineInputs, {
    locale?: "en" | "de";
}, {}>;
