/* eslint-disable */
import * as registry from '../registry.js'
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ count: NonNullable<unknown> }} Settings_Filter_CountInputs */

const en_settings_filter_count = /** @type {(inputs: Settings_Filter_CountInputs) => LocalizedString} */ (i) => {const countPlural = registry.plural("en", i?.count, {}, true);
	if (countPlural === "one") return /** @type {LocalizedString} */ (`${i?.count} setting matches`);
	if (countPlural === "other") return /** @type {LocalizedString} */ (`${i?.count} settings match`);
	return /** @type {LocalizedString} */ ("settings_filter_count");
};

const de_settings_filter_count = /** @type {(inputs: Settings_Filter_CountInputs) => LocalizedString} */ (i) => {const countPlural = registry.plural("de", i?.count, {}, true);
	if (countPlural === "one") return /** @type {LocalizedString} */ (`${i?.count} Einstellung passt`);
	if (countPlural === "other") return /** @type {LocalizedString} */ (`${i?.count} Einstellungen passen`);
	return /** @type {LocalizedString} */ ("settings_filter_count");
};

/**
* | countPlural | output |
* | --- | --- |
* | "one" | "{count} setting matches" |
* | "other" | "{count} settings match" |
*
* @param {Settings_Filter_CountInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const settings_filter_count = /** @type {((inputs: Settings_Filter_CountInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_CountInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_settings_filter_count(inputs)
	return en_settings_filter_count(inputs)
});