/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Search_EmptyInputs */



/**
* | output |
* | --- |
* | "No pages match the search." |
*
* @param {Docs_Search_EmptyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_search_empty = /** @type {((inputs?: Docs_Search_EmptyInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_EmptyInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Keine Seite passt zur Suche.")
	return /** @type {LocalizedString} */ ("No pages match the search.")
});