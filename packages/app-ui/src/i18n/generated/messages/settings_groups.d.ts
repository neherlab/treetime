export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_GroupsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_GroupsInputs */
/**
* | output |
* | --- |
* | "Setting groups" |
*
* @param {Settings_GroupsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_groups: ((inputs?: Settings_GroupsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_GroupsInputs, {
    locale?: "en" | "de";
}, {}>;
