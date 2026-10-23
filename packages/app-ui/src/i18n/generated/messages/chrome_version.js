/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ version: NonNullable<unknown> }} Chrome_VersionInputs */

const en_chrome_version = /** @type {(inputs: Chrome_VersionInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`v${i?.version}`)
};

const de_chrome_version = /** @type {(inputs: Chrome_VersionInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`v${i?.version}`)
};

/**
* | output |
* | --- |
* | "v{version}" |
*
* @param {Chrome_VersionInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_version = /** @type {((inputs: Chrome_VersionInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_VersionInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_chrome_version(inputs)
	return en_chrome_version(inputs)
});