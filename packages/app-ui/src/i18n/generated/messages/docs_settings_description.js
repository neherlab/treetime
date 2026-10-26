/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ command: NonNullable<unknown>, title: NonNullable<unknown> }} Docs_Settings_DescriptionInputs */

const en_docs_settings_description = /** @type {(inputs: Docs_Settings_DescriptionInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Settings of treetime ${i?.command} (${i?.title}), with their flags, config keys, and defaults.`)
};

const de_docs_settings_description = /** @type {(inputs: Docs_Settings_DescriptionInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Einstellungen von treetime ${i?.command} (${i?.title}), mit ihren Optionen, Konfigurationsschlüsseln und Standardwerten.`)
};

/**
* | output |
* | --- |
* | "Settings of treetime {command} ({title}), with their flags, config keys, and defaults." |
*
* @param {Docs_Settings_DescriptionInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_settings_description = /** @type {((inputs: Docs_Settings_DescriptionInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Settings_DescriptionInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_docs_settings_description(inputs)
	return en_docs_settings_description(inputs)
});