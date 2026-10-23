export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_ConflictsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Settings_ConflictsInputs */
/**
* | output |
* | --- |
* | "Cannot be combined with" |
*
* @param {Settings_ConflictsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_conflicts: ((inputs?: Settings_ConflictsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_ConflictsInputs, {
    locale?: "en" | "de";
}, {}>;
