/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ value: NonNullable<unknown> }} Settings_DefaultInputs */

const en_settings_default = /** @type {(inputs: Settings_DefaultInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Default: ${i?.value}`)
};

const de_settings_default = /** @type {(inputs: Settings_DefaultInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Standard: ${i?.value}`)
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
export const settings_default = /** @type {((inputs: Settings_DefaultInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_DefaultInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_settings_default(inputs)
	return en_settings_default(inputs)
});