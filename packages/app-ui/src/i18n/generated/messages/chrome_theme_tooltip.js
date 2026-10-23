/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{ label: NonNullable<unknown> }} Chrome_Theme_TooltipInputs */

const en_chrome_theme_tooltip = /** @type {(inputs: Chrome_Theme_TooltipInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`${i?.label}. Click to change.`)
};

const de_chrome_theme_tooltip = /** @type {(inputs: Chrome_Theme_TooltipInputs) => LocalizedString} */ (i) => {
	return /** @type {LocalizedString} */ (`${i?.label}. Zum Wechseln klicken.`)
};

/**
* | output |
* | --- |
* | "{label}. Click to change." |
*
* @param {Chrome_Theme_TooltipInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_theme_tooltip = /** @type {((inputs: Chrome_Theme_TooltipInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_TooltipInputs, { locale?: "en" | "de" }, {}>} */ ((inputs, options = {}) => {
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return de_chrome_theme_tooltip(inputs)
	return en_chrome_theme_tooltip(inputs)
});