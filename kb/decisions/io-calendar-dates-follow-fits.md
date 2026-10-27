# Calendar dates of numeric dates follow FITS and ISO 8601

TreeTime writes the calendar date of every numeric date inside the range of chrono. Dates before the year 1 and after the year 9999 use the signed five-digit year of the <a id="gloss-use-1"></a>FITS <sup>[1](#gloss-1)</sup> standard, a fixed form of the ISO 8601 expanded year. The month and day of a date outside the years 1 to 9999 can differ from TreeTime 0.x by up to one day.

**Type**: Scientific output change from v0.

**Status**: Approved by the maintainer team.

## Numeric date

A numeric date is the year plus the part of that year that passed since 1 January 00:00, with the real length of the year (365 or 366 days):

$$
t = y + \frac{s}{S_y}
$$

where:

- $t$ -- numeric date
- $y$ -- calendar year, $y = \lfloor t \rfloor$
- $s$ -- seconds since 1 January 00:00 of year $y$
- $S_y$ -- seconds in year $y$

A calendar day maps to its noon, $t = y + (d - 0.5) / D_y$, where $d$ is the day of the year and $D_y$ the days in year $y$. This is the definition of v0 `numeric_date()` [`packages/legacy/treetime/treetime/utils.py#L146-L164`](../../packages/legacy/treetime/treetime/utils.py#L146-L164) and of the `decimalyear` format of astropy [[src](https://github.com/astropy/astropy/blob/bf1d238b4da64df38eb69695f0d36275fe0d6292/astropy/time/formats.py#L658)]. v1 takes the year with $\lfloor t \rfloor$ also for negative values, so the part of the year is never negative.

## Calendar

The calendar is the proleptic Gregorian calendar of chrono with astronomical year numbering: year 0 is 1 BCE, year -1 is 2 BCE. chrono covers the years -262143 to 262142 (`NaiveDate::MIN` and `NaiveDate::MAX`). All leap-year and month-length facts come from chrono. A value outside this range, NaN, and the infinities have no calendar date.

## String form

- **Years 0 to 9999**: `YYYY-MM-DD`, for example `2012-06-15` and `0000-07-02`
- **Other years**: the year with a mandatory sign and at least five digits, `{year:+06}-MM-DD`, for example `-01910-05-30`, `-00005-01-01`, and `+10000-07-02`. Years beyond ±99999 get six digits, for example `-262143-07-02`
- **No calendar date**: the dates table writes `--`, the augur node data has `numdate` and no `date`, and the app shows the numeric date with two decimals
- **No warning**: a root before the year 1 gives no run warning

The FITS Standard 4.0 <a id="cite-1"></a>[IAU FITS Working Group 2016](https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf) [[1](#ref-1)] defines the date strings of FITS in its section 9.1.1 as a subset of ISO 8601: the unsigned four-digit year `CCYY-MM-DD`, or the signed five-digit year `[±C]CCYY-MM-DD`. The same section counts years with a year zero (year 0 = 1 BCE), reads dates before 1582 in the proleptic Gregorian calendar, and limits the five-digit form to the years -99999 to +99999. The signed five-digit year came into FITS with the time coordinates of <a id="cite-2"></a>[Rots et al. 2015](https://doi.org/10.1051/0004-6361/201424653) [[2](#ref-2)]. ISO 8601 allows a year with more digits and a sign when the reader and the writer agree on the number of digits; FITS sets this number to five. The six-digit years beyond ±99999 are outside FITS and inside the expanded form of ISO 8601.

astropy implements this standard as the `longdate` subformat of its `fits` time format [[src](https://github.com/astropy/astropy/blob/bf1d238b4da64df38eb69695f0d36275fe0d6292/astropy/time/formats.py#L2172-L2213)]. It reads the form back only for the years -4799 to 99999: its pattern accepts exactly five year digits, and ERFA rejects years before -4799 [[src](https://github.com/liberfa/erfa/blob/1a8044cde5b7763295d472a6443387239127c6c8/src/cal2jd.c#L57)].

## Rationale

- Dates before the year 1 occur on deep trees, on ancient samples, and when a weak clock signal puts the root far into the past
- One-year accuracy is enough for these dates; the calendar date keeps the column readable, and the numeric date stays the value for computation
- The `YYYY-MM-DD` shape of the years 0 to 9999 stays unchanged, and the FITS form is a published standard that tools can parse

## Difference from v0

v0 `datestring_from_numeric()` [`packages/legacy/treetime/treetime/utils.py#L196-L214`](../../packages/legacy/treetime/treetime/utils.py#L196-L214) writes a date outside the years 1 to 9999 with the year $\lfloor t \rfloor$, the month and day of the same part of the year in 1900 (a common year), and the format `'%04d-%02d-%02d'`.

- **Year form**: v0 writes `-1910-05-30` for $t = -1909.59$ and `-005-...` for the year -5. v1 writes `-01910-05-30` and `-00005-...`
- **Month and day**: v1 takes them from the real year, v0 from 1900. In a leap year outside the years 1 to 9999, the month and day can differ from v0 by one day
- **Years 1 to 9999**: no change from this decision. The one-day difference of v0 at the boundary of a day is [kb/issues/N-io-node-data-date-string-fp-boundary.md](../issues/N-io-node-data-date-string-fp-boundary.md)

## Numeric input dates

A numeric date in the metadata keeps its value as written, as in v0, also outside the years 1000 to 3000 (`-500` stays `-500.0`). A compact calendar date such as `20200315` still reads as 15 March 2020, because the date formats are tried before the number.

## Implementation

[`packages/treetime-utils/src/datetime/year_fraction.rs`](../../packages/treetime-utils/src/datetime/year_fraction.rs): `fn year_fraction_to_date()`, `fn year_fraction_to_datestring()`, and `fn year_fraction_days_between()` return `None` for a value without a calendar date.

## Glossary

1. <a id="gloss-1"></a> **FITS (Flexible Image Transport System).** The data format standard of astronomy, maintained by the IAU FITS Working Group. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> IAU FITS Working Group. 2016. _Definition of the Flexible Image Transport System (FITS): The FITS Standard, Version 4.0_. International Astronomical Union. https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf [↩](#cite-1)
2. <a id="ref-2"></a> Rots, Arnold H., Peter S. Bunclark, Mark R. Calabretta, Steven L. Allen, Richard N. Manchester, and William T. Thompson. 2015. "Representations of Time Coordinates in FITS." _Astronomy & Astrophysics_ 574:A36. https://doi.org/10.1051/0004-6361/201424653 [↩](#cite-2)
