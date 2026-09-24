[![GitHub Workflow Status (building)](https://img.shields.io/github/actions/workflow/status/soilwater/crnpy/python-package.yml)](https://github.com/soilwater/crnpy/actions/workflows/python-package.yml)
[![GitHub Workflow Status (publish)](https://img.shields.io/github/actions/workflow/status/soilwater/crnpy/python-publish.yml?label=publish)](https://github.com/soilwater/crnpy/actions/workflows/python-publish.yml)
[![PyPI - Status](https://img.shields.io/pypi/v/crnpy)](https://pypi.org/project/crnpy/)
[![GitHub commits since latest release (by SemVer including pre-releases)](https://img.shields.io/github/commits-since/soilwater/crnpy/latest/main)](https://github.com/soilwater/crnpy)
[![JOSS submission status](https://joss.theoj.org/papers/e65c1bb5fee58c39289efc4547d1fd10/status.svg)](https://joss.theoj.org/papers/e65c1bb5fee58c39289efc4547d1fd10)

# Cosmic-Ray Neutron Python (CRNPy) Library

<img src="https://raw.githubusercontent.com/soilwater/crnpy/main/docs/img/logo/crnpy-logo.png" alt="CRNPY logo" width="250"/>

## Overview
Welcome to the homepage of the CRNPy (Cosmic-Ray Neutron Python) library, an open-source Python library designed for the processing and conversion of raw neutron counts from cosmic-ray neutron probes (CRNP) into soil moisture data.

This library has been developed with the intent of providing a comprehensive yet easy-to-use workflow for processing raw data from a variety of CRNP, encompassing multiple manufacturers and models.

## Statement of Need
CRNPs are a valuable tool for non-invasive soil moisture estimation at the hectometer scale (e.g., typical agricultural fields), filling the gap between point-level sensors and large-scale (i.e., several kilometers) remote sensors onboard orbiting satellites. However, cleaning, processing, and analyzing CRNP data involves multiple corrections and filtering steps spread across multiple peer-reviewed manuscripts. CRNPy simplifies these steps by providing a complete, user-friendly, and well-documented library with minimal dependencies that includes examples to convert raw CRNP data into soil moisture. The library is designed to be accessible to both researchers and instrument manufacturers. Unlike other similar libraries, CRNPy does not require any specific naming convention for the input data or large external data sources, or reanalysis data.

## Key Features
- Versatile and instrument agnostic: CRNPy can handle data from various CRNP manufacturers and models. It has been successfully tested on both roving and stationary CRNP.

- Modular: The library is designed to be modular, allowing users to easily customize the processing workflow to their needs.


## Installation
To install the CRNPy library, you can use Python's package manager. Open a terminal and type:

```pip install crnpy```

from the Jupyter notebook, type:

```!pip install crnpy```

Ideally dependencies should be installed automatically. If not, you can install them manually by typing:

```pip install -r requirements.txt```

The CRNPy library is compatible with Python 3.8 and above.
See [requirements.txt](https://github.com/soilwater/crnpy/blob/main/requirements.txt) for a list of dependencies.

## Examples

- [https://soilwater.github.io/crnpy/examples/stationary/example_RDT_station/](Processing and analyzing data from a stationary detector)
- [https://soilwater.github.io/crnpy/examples/rover/Hydroinnova_rover_example/](Processing and analyzing data from a roving detector)
- [https://soilwater.github.io/crnpy/examples/calibration/calibration/](Device-specific field calibration)

## Authors
The CRNPy library was developed at the Kansas State University Soil Water Processes Lab by:

- Joaquin Peraza

- Andres Patrignani

The Soil Water Processes Lab at Kansas State University combines a range of experimental and computational approaches to tackle pressing issues in soil and water research. The development of the CRNPy library is a step forward to creating reproducible data processing workflows across the scientific community using cosmic-ray neutrons probes for soil moisture sensing.


## Community Guidelines
### Contributing
To contribute to the software, please first fork the repository and create your own branch from `main`. Ensure your code adheres to our established code structure and includes appropriate test/examples coverage. CRNPy source code is located in the `/src/crnpy/` folder, and tests implemented using `pytest` are stored in the `/src/tests/` folder. Submit a pull request with a clear and detailed description of your changes to include them in the main repository.

To build the documentation locally, install the documentation toolchain with `pip install -r requirements-docs.txt` and run `mkdocs build` (or `mkdocs serve`) from the repository root.

### Reporting Issues
If you encounter any issues or problems with the software, please report them on our [issues page](https://github.com/soilwater/crnpy/issues). Include a detailed description of the issue, steps to reproduce the problem, any error messages you received, and details about your operating system and software version.

### Seeking Support
If you need support, please first refer to the documentation. If you still require assistance, post a question on the [issues page](https://github.com/soilwater/crnpy/issues) with the `question` tag. For private inquiries, you can reach us via email at jperaza@ksu.edu or andrespatrignani@ksu.edu.

## Changelog

### Version 0.8.0
Changes since version 0.7.0. Items marked **breaking** change the public API.

- **breaking** `get_reference_neutron_flux()` now requires an explicit `date` argument (previously defaulted to 2011-05-01). This prevents silently requesting a reference date at which the selected station may have no data (e.g. stations that came online after 2011-05-01). A value of 2011-05-01 remains the conventional choice, following Zreda et al. (2012), Hawdon et al. (2014) and Bogena et al. (2022).
- `get_incoming_neutron_flux()` now distinguishes NMDB's definitive "no data available" response, which is reported once and returned as None without retrying, from transient failures (empty or unexpected pages), which are retried up to three times. A 30 s request timeout and handling of network errors were added, and a new `report_no_data` argument controls the no-data message.
- `find_neutron_monitor()`, when given a period, now lists only the monitors that have data for that period (the tentative stations without data are skipped silently, since the user did not choose them), returns a DataFrame with a reset index and without the redundant `Period available` column, and prints the results once instead of both printing and returning them.
- `correction_road()` now exposes all eleven parameters of Schrön et al. (2018): the second weight of the distance term is an independent `p9` (default 0.06) and its decay rate is `p10` (default 0.01). Results are unchanged at the default parameters.
- `interpolate_incoming_flux()` gains a `tolerance` argument (default 1 hour) so CRNP timestamps with no NMDB observation within the tolerance remain NaN instead of taking a far-away value.
- `sensing_depth()` raises a clear `ValueError` when `method='Schron_2017'` and `dist` is not provided, instead of failing with an obscure error.
- `correction_incoming_flux()`: `fill_na` now works for NumPy array input in addition to pandas Series, and the default reference flux is taken from the first value by position, so a pandas Series with a non-default index no longer raises an error.
- `exp_filter()` returns the input unchanged for a constant series instead of dividing by zero.
- `idw()` (and `interpolate_2d(method='idw')`) no longer emit a `RuntimeWarning` for grid cells with no observation within `neighborhood`; such cells are left as NaN, as before, but without the spurious 0/0 division.

### Version 0.7.0
Changes since version 0.6.1. Items marked **breaking** change the public API.

- **breaking** `nrad_weight()` now implements only the revised footprint weighting of Schrön et al. (2017). The Köhli et al. (2015) weighting and the `method` argument were removed; `profiles` and `p` must be provided. Sample `depth` is in cm, consistent with the penetration depth D86, which the paper defines in cm.
- Fixed `nrad_weight()`: the horizontal weighting function was evaluated with soil moisture in place of air humidity and vegetation height in place of soil moisture, so air humidity never entered the weights and profiles beyond 50 m received negative weights. Both weighting functions now use the site air humidity and the field-average soil moisture of the current iteration, as in Sect. 3 of Schrön et al. (2017) and its R/MATLAB supplement. Field-average values from earlier versions will differ.
- Fixed `correction_incoming_flux()` with `Rc_method='McJannetandDesilets2023'`: the factor returned was the reciprocal of what the library's convention (counts divided by `fi`) requires, so the correction acted in the wrong direction. It now returns `tau*(I/Iref) + 1 - tau` following Eq. 10 of McJannet and Desilets (2023), and reduces to `I/Iref` when site and reference coincide.
- Fixed `correction_road()`: the moisture term had the wrong sign on `p2` and omitted `p5`, which inflated instead of reduced the counts. It now implements Eq. 6 and Table 1 of Schrön et al. (2018) with a new `p5=0.39` argument, and a road width of zero returns the counts unchanged instead of raising an error.
- `correction_bwe()`: the default `r2_N0` was 0.05, ten times the value of Baatz et al. (2015). The default is now 0.0053 (r2 = 6.4 cph per kg m-2 BWE, N0 = 1210 cph).
- Fixed `uncertainty_counts(metric='cv')`: the coefficient of variation was multiplied by the correction factors; it is now `1/sqrt(N)` (Jakobi et al., 2020, Eq. 7).
- `sensing_depth()`: the documentation now states the result is in cm for both methods (Franz et al., 2012, Eq. 5; Schrön et al., 2017, Eq. 4). In the `Schron_2017` method the lattice water is converted to volumetric water equivalent with the bulk density before entering D86 (Schrön et al., 2017, Eq. 2), as already done in the `Franz_2012` method; this changes the Schrön depth by a few percent.
- `cutoff_rigidity()`: the bundled world grid was replaced by the published grid of calculated vertical cutoff rigidities for epoch 1995.0 of Smart and Shea (2008, Proc. 30th ICRC, 1, 733-736), tabulated every 5 degrees in latitude and 15 degrees in longitude, and west longitudes are now mapped onto the 0-360 east longitude grid (previously they were mirrored, which happened to compensate a shift in the western half of the old table). Against 102 neutron monitors the mean absolute error drops from 0.27 to 0.18 GV. Returned values change by up to a few tenths of a GV, and by more in South America. The function now cites only that source.
- `lattice_water()`: the documentation now states that inputs and output are gravimetric (percent by mass), and how to convert the result to the g/g fraction used as `Wlat` elsewhere in the library. No change in the calculation.
- Fixed `exp_filter()`: a single missing value made the rest of the series NaN. The gain now uses the time elapsed since the last available observation (Albergel et al., 2008, Eq. 6), so gaps are skipped and the filter continues. Results are unchanged for series without gaps.
- Fixed `remove_incomplete_intervals()`: the first row was always removed because its time difference is undefined; it is now kept unless `remove_first=True`.
- Fixed `smooth_1d()` with a DataFrame and `method='savitzky_golay'`, which raised an error and modified the input in place.
- Fixed `is_outlier()` when `min_val`/`max_val` are omitted (it raised a type error). The `'range'` method listed in the documentation is now implemented, and `'scaled_mad'` is two-sided (values below the median are also flagged), as in MATLAB `isoutlier`.
- Fixed `find_neutron_monitor()` when no station has data for the requested period (it raised an error); the ten closest stations are now returned with `Period available` set to False.
- `get_incoming_neutron_flux()` tries at most three times, pausing five seconds between attempts, when the NMDB endpoint returns a page without the data block, which it does intermittently even when data are available. If no data block is found after three attempts it prints a message suggesting to wait a minute and returns None cleanly (previously a single failed response returned None, which made `find_neutron_monitor()` report every station as unavailable during such spells).
- `interpolate_incoming_flux()` keeps periods without NMDB data as NaN instead of failing an assertion; use `fill_na` in `correction_incoming_flux()`.
- `idw()` returns the observed value at prediction points that coincide with an observation instead of NaN.
- `total_raw_counts()` fills a missing detector with the mean of the other detectors only when more than one detector column is present (the previous check looked at the number of rows).
- Error messages raised as plain strings in `spatial_average()`, `interpolate_2d()` and `uncertainty_counts()` are now proper `ValueError` exceptions; the NaN warning in `interpolate_2d()` now reports the actual counts.
- Compatibility with pandas 2.2 and 3: `fill_missing_timestamps()` defaults to the frequency alias `'h'` (the `'H'` alias is deprecated in pandas 2.2 and removed in pandas 3) and the examples use it; datetime columns of any resolution or time zone are accepted; chained in-place assignments and column assignments on filtered frames, which warn in pandas 2 and stop working in pandas 3, were removed from the examples and tests.
- All references in the function docstrings were standardized to APA 7th style (authors as `Last, F. M.` with `&` before the final author, year after the authors, sentence-case titles, full journal names, en-dash page ranges, DOIs as `https://doi.org/…`). Existing DOIs were kept; none were added to references that lacked one.
- Added `requirements-docs.txt` with the pinned documentation toolchain (MkDocs 1.x, Material, mkdocstrings, mkdocs-autorefs, mkdocs-jupyter, glightbox); the documentation page on correction routines had six cross-references to non-existent function names, now fixed, and describes the footprint weighting used for calibration.
- Documentation corrections in `correction_pressure()` (reference pressure and attenuation length), `correction_humidity()` (input is absolute humidity), `counts_to_vwc()` (full equation and gravimetric units of `Wlat` and `Wsoc`, after Hawdon et al., 2014, Eq. 7), `abs_humidity()` (returns absolute humidity), `get_incoming_neutron_flux()` (NMDB acknowledgement, printed with `verbose=True`) and `latlon_to_utm()`.

