---
title: 'CRNPY: An Open-Source Python Library for Cosmic Ray Neutron Probe Data Processing'
tags:
  - Python
  - Cosmic Ray Neutron Probes
  - Soil Moisture
  - Data Processing
authors:
  - name: Joaquin Peraza
    orcid: 0000-0002-XXXX-XXXX
    affiliation: 1
  - name: Andres Patrignani
    orcid: 0000-0002-XXXX-XXXX
    affiliation: 1
affiliations:
 - name: Kansas State University, Department of Agronomy
   index: 1
date: 01 June 2023
bibliography: paper.bib

# Summary

CRNPY is a Python library developed for processing and analyzing cosmic ray neutron
probe (CRNP) data. The library provides a comprehensive set of tools for handling raw
CRNP data, including atmospheric and calibration corrections, data filtering, and soil
moisture estimation. The library is designed to be user-friendly and flexible, 
accommodating a wide range of use cases and data formats.


# Statement of Need

Cosmic ray neutron probes are a valuable tool for estimating soil moisture at the field
scale filling the gap between accurate point level estimates and medium to high scale
remote sensing and modelling products. However, processing and analyzing CRNP data can
be complex and time-consuming. CRNPY simplifies this process, providing a user-friendly
and efficient way to transform raw CRNP data into actionable information.
The library is designed to be accessible to both researchers and practitioners, 
with clear documentation and examples provided.

# Methodology

CRNPY implements a range of methods for processing and analyzing CRNP data. These include:

- Atmospheric and calibration corrections: CRNPY provides functions for correcting raw neutron counts for atmospheric pressure, humidity, and temperature variations, as well as for incoming neutron flux [@Zreda2012; @Andreasen2017].
- Data filtering: The library includes functions for removing outliers and periods of invalid data from the raw neutron counts.
- Soil moisture estimation: CRNPY implements the Desilets method for estimating soil moisture from corrected neutron counts [@Desilets2010].
- Biomass correction: The library provides a function for correcting neutron counts for the effects of above-ground biomass, following the approach described in Baatz et al. (2015).

# References

The library is based on the methods described in the following references:

- Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., et al. (2012). COSMOS: the cosmic-ray soil moisture observing system. Hydrol. Earth Syst. Sci. 16, 4079–4099
- Andreasen, M., Jensen, K. H., Desilets, D., Franz, T. E., Zreda, M., Bogena, H. R., et al. (2017). Status and perspectives of cosmic-ray neutron sensors for soil moisture monitoring and neutron detector development. Vadose Zone J. 16.
- Desilets, D., Zreda, M., & Ferré, T. P. A. (2010). Nature's neutron probe: Land surface hydrology at an elusive scale with cosmic rays. Water Resources Research, 46(11).
- Baatz, R., Bogena, H. R., Hendricks Franssen, H. J., Huisman, J. A., Qu, W., Montzka, C., et al. (2015). An empirical vegetation correction forI apologize for the oversight. I've removed the installation section and made some adjustments to better align with the JOSS submission guidelines. Here's the revised version:
