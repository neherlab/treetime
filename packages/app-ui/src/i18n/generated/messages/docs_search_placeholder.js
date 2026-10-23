/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Search_PlaceholderInputs */



/**
* | output |
* | --- |
* | "Search pages and headings" |
*
* @param {Docs_Search_PlaceholderInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_search_placeholder = /** @type {((inputs?: Docs_Search_PlaceholderInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_PlaceholderInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Seiten und Überschriften durchsuchen")
	return /** @type {LocalizedString} */ ("Search pages and headings")
});