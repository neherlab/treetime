/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ current: NonNullable<unknown>, preferred: NonNullable<unknown> }} Docs_Version_BannerInputs */

const en_docs_version_banner = /** @type {(inputs: Docs_Version_BannerInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`These are the docs of TreeTime ${i?.current}. The recommended version is ${i?.preferred}.`)
};

const de_docs_version_banner = /** @type {(inputs: Docs_Version_BannerInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Dies ist die Doku von TreeTime ${i?.current}. Empfohlen ist die Version ${i?.preferred}.`)
};

/**
* | output |
* | --- |
* | "These are the docs of TreeTime {current}. The recommended version is {preferred}." |
*
* @param {Docs_Version_BannerInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_version_banner = /** @type {((inputs: Docs_Version_BannerInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Version_BannerInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_docs_version_banner(inputs)
	return en_docs_version_banner(inputs)
});