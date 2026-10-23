/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Search_ResultsInputs */



/**
* | output |
* | --- |
* | "Search results" |
*
* @param {Docs_Search_ResultsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_search_results = /** @type {((inputs?: Docs_Search_ResultsInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_ResultsInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Suchergebnisse")
	return /** @type {LocalizedString} */ ("Search results")
});