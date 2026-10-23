/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Chrome_GithubInputs */



/**
* | output |
* | --- |
* | "TreeTime on GitHub" |
*
* @param {Chrome_GithubInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_github = /** @type {((inputs?: Chrome_GithubInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_GithubInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("TreeTime auf GitHub")
	return /** @type {LocalizedString} */ ("TreeTime on GitHub")
});