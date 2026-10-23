/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ preferred: NonNullable<unknown> }} Docs_Version_Banner_LinkInputs */

const en_docs_version_banner_link = /** @type {(inputs: Docs_Version_Banner_LinkInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Open this page in ${i?.preferred}`)
};

const de_docs_version_banner_link = /** @type {(inputs: Docs_Version_Banner_LinkInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`Diese Seite in ${i?.preferred} öffnen`)
};

/**
* | output |
* | --- |
* | "Open this page in {preferred}" |
*
* @param {Docs_Version_Banner_LinkInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_version_banner_link = /** @type {((inputs: Docs_Version_Banner_LinkInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Version_Banner_LinkInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_docs_version_banner_link(inputs)
	return en_docs_version_banner_link(inputs)
});