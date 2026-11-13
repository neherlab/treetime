export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Command_Line_Left_OutInputs = {
    flag: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "{flag} is left out for this value" |
*
* @param {Settings_Command_Line_Left_OutInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_command_line_left_out: ((inputs: Settings_Command_Line_Left_OutInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Command_Line_Left_OutInputs, {
    locale?: "en" | "de";
}, {}>;
