---
title: "CRNPy: An Open-Source Python Library for Cosmic Ray Neutron Probe Data Processing"
tags:
- Python
- Cosmic Ray Neutron Probes
- Soil Moisture
- Data Processing
date: "01 June 2023"
output:
  pdf_document:
    fig_caption: yes
  html_document:
    df_print: paged
  word_document: default
authors:
- name: Peraza Rud, Joaquin
  orcid: "0009-0007-9337-830X"
  affiliation: 1
- name: Tyson Ochsner
  orcid: "0000-0003-0875-4491"
  affiliation: 2
- name: Andres Patrignani
  orcid: "0000-0001-5955-5877"
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Department of Agronomy, Kansas State University, Manhattan, KS, USA.
  index: 1
- name: Department of Plant and Soil Sciences, Oklahoma State University, Stillwater, OK, USA.
  index: 2
editor_options:
  markdown:
    wrap: sentence
---

# Summary

CRNPy is a Python library developed for processing and analyzing cosmic-ray neutron counts recorded by stationary and roving moderated neutron probes (CRNP) for soil moisture sensing. The CRNPy library includes common routines for implementing atmospheric corrections, biomass corrections, road corrections, one- and two-dimensional filtering, horizontal and vertical footprint analysis, uncertainty estimation, depth extrapolation operators, device calibration, and conversion of neutron counts into volumetric soil water content. The library is designed to be user-friendly, versatile, and instrument agnostic to enable users the integration of the CRNPy library in a wide range of agricultural and hydrological applications.

# Statement of Need

Cosmic ray neutron probes are a valuable tool for non-invasive soil moisture sensing at the hectometer scale (e.g., typical agricultural fields), filling the gap between point-level sensors and large-scale (i.e., several kilometers) remote sensors onboard orbiting satellites. However, cleaning, processing, and analyzing CRNP data involves multiple corrections and filtering steps spread across multiple peer-reviewed articles. CRNPy simplifies these steps by providing a complete, user-friendly, and well-documented library with minimal dependencies that includes real-world examples to convert raw CRNP data into soil moisture. The library is designed to be accessible to both researchers and instrument manufacturers seeking to integrate the library within their hardware. Unlike other similar libraries, CRNPy does not require any specific naming convention for the input data or the download of large external data sources or reanalysis data.

# Library features

The CRNPy library includes common routines for creating workflows to process raw data from both stationary (\autoref{fig:workflow_stationary}) and roving (\autoref{fig:workflow_rover}) devices. Some of the notable features of the CRNPy library include:

- The entire CRNPy library is implemented on Python programming language using, Numpy, Pandas, SciPy, and Matplotlib libraries, which are included in common Python data science bundles like the Anaconda open source ecosystem.
- Utility functions for determining site-specific information required before processing raw neutron counts, including the determination of lattice water (i.e., bounded water to clay particles), geomagnetic cutoff rigidity [@smart2001geomagnetic], and searching for a reference neutron monitor [@2009NMDB].
- Delimited-text files can be used without the need to change column names, which increases reproducibility and minimizes human error. Each function of the CRNPy library accepts either a Numpy array or a Pandas Series, enabling a more versatile, modular, and customizable workflow that adapts to instrument outputs from different manufacturers. The CRNPy examples using multiple devices provide guidelines and a starting point for users getting started with CRNP devices.
- Outliers can be caused by hardware malfunction, vibrations, or external neutron sources (e.g., use of nearby neutron probe soil moisture meters like the 503 Hydroprobe by Instrotek, Inc.). The CRNPy library offers several methods for outliers detection using the `is_outlier()` method, which includes a simple range detection based on user-provided lower and upper boundaries, interquartile range, z-scores, and a scaled mean absolute difference [@iglewicz1993volume].
- The ability to compute the total raw neutron counts for CRNP devices with multiple neutron detectors.
- Implement corrections for atmospheric pressure, air humidity [@zreda2012cosmos; @andreasen2017status], and incoming neutron flux [@rosolem2013] using the `correction_pressure()`, `correction_humidity()`, and `correction_incoming_flux()` functions. This step is necessary to account for the effects of atmospheric conditions on neutron attenuation. \autoref{fig:output_stationary}a and \autoref{fig:output_stationary}b show the impact of each correction in the raw neutron count observed with a stationary CRNP.
- Additional corrections may be necessary to account for the effects of other surrounding hydrogen pools, such as above- and below-ground plant biomass and road conditions, on the neutron counts. `biomass_to_bwe()` and `correction_bwe()` functions are used for biomass correction [@baatz2015empirical; @wahbi2018situ], while the road_correction() function is used for correcting for differences in road moisture and surrounding field soil moisture when doing mobile surveys [@schron2018cosmic]. These steps are optional and depend on the specific research needs and data conditions.
- The conversion of corrected counts into volumetric soil water content follows the approach suggested by [@desilets2010nature] using the `counts_to_vwc()` function.
- Neutron count uncertainty can estimated with the `uncertainty_counts()` function and propagated into the resulting volumetric water content using the `uncertainty_vwc()` function by following the method detailed in @jakobi2020error. \autoref{fig:output_stationary}c shows the neutron count uncertainty propagated to the volumetric water content.
- The `sensing_depth()` method allows the sensing depth estimation of the measured soil volumetric water content, following the approach of previous studies finding the volume that accounts for 86% (2 e-folds assuming an exponential decay) of the origin of the counted neutrons [@zreda2008measuring; @franz2012measurement; @kohli2015footprint].
- Because observations with CRNP devices typically represent the soil moisture conditions in the top 10-20 cm, an exponential filter operator [@albergel2008near] is provided to extrapolate soil moisture conditions at greater soil depths using the function `exp_filter()`, which can be used to extend soil moisture conditions in the rootzone [@rossini2021predicting; @franz2020practical]. \autoref{fig:output_stationary}d shows the variation of the sensing depth over time for a stationary CRNP.
- For roving devices, CRNPy includes a few utility functions for spatial filtering and basic interpolation routines with cubic, linear, nearest neighbor, and inverse distance weighting interpolation methods. \autoref{fig:output_rover} represents the output of CRNP rover transect, after being processed and interpolated using CRNPy.

![Example workflow for stationary CRNP, dashed lines represent optional steps.\label{fig:workflow_stationary}](figures/workflow_rdt.png)

![Example workflow for roving CRNP, dashed lines represent optional steps.\label{fig:workflow_rover}](figures/workflow_hydroinnova.png)

![Outputs of a stationary device in each of the steps of the workflow, A) Impact of each factor in the correction process. B) Difference between raw and corrected counts. C) Resulting volumetric water content time series with the propageted neutron count uncertainity. D) Depth that accounts for the 86% of the observed counts.  \label{fig:output_stationary}](figures/timeseries.png)

![Contour plot of the roving transect spatial output. Each marker represents the estimated center of the observation.\label{fig:output_rover}](figures/rover.png)

# Acknowledgements

This project was supported by Agriculture and Food Research Initiative Competitive Grant no. 2019-68012-29888 from the USDA National Institute of Food and Agriculture. The authors declare no competing interests.

# References

