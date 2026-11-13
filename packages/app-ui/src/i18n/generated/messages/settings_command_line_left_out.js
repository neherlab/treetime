/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ flag: NonNullable<unknown> }} Settings_Command_Line_Left_OutInputs */

const en_settings_command_line_left_out = /** @type {(inputs: Settings_Command_Line_Left_OutInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`${i?.flag} is left out for this value`)
};

const de_settings_command_line_left_out = /** @type {(inputs: Settings_Command_Line_Left_OutInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`${i?.flag} entfällt bei diesem Wert`)
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
export const settings_command_line_left_out = /** @type {((inputs: Settings_Command_Line_Left_OutInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Command_Line_Left_OutInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_settings_command_line_left_out(inputs)
	return en_settings_command_line_left_out(inputs)
});