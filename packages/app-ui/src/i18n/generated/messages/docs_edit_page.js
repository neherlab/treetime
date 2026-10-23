/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Edit_PageInputs */



/**
* | output |
* | --- |
* | "Edit this page on GitHub" |
*
* @param {Docs_Edit_PageInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_edit_page = /** @type {((inputs?: Docs_Edit_PageInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Edit_PageInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Diese Seite auf GitHub bearbeiten")
	return /** @type {LocalizedString} */ ("Edit this page on GitHub")
});