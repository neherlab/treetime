/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ label: NonNullable<unknown> }} Settings_AboutInputs */

const en_settings_about = /** @type {(inputs: Settings_AboutInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`About ${i?.label}`)
};

const de_settings_about = /** @type {(inputs: Settings_AboutInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Über ${i?.label}`)
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
export const settings_about = /** @type {((inputs: Settings_AboutInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_AboutInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_settings_about(inputs)
	return en_settings_about(inputs)
});