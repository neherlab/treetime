/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_SearchInputs */



/**
* | output |
* | --- |
* | "Search the docs" |
*
* @param {Docs_SearchInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_search = /** @type {((inputs?: Docs_SearchInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_SearchInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Doku durchsuchen")
	return /** @type {LocalizedString} */ ("Search the docs")
});