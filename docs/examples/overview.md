# CRNPy Examples Overview

This document provides an overview of the examples available for the Cosmic Ray Neutron Python (CRNPy) library. The examples demonstrate how to process and analyze neutron detection data from various detectors.

## [Stationary Example: RDT Station](../stationary/example_RDT_station/)

This example demonstrates how to process and analyze neutron detection data from 3 Radiation Detection Technologies, Inc (RDT) neutron detectors using the CRNPy library. The tutorial covers steps including data loading, count data processing, normalization, outlier removal, atmospheric correction, and conversion of neutron counts to Volumetric Water Content (VWC).

Key steps include:

- Loading data from a .csv file and converting the 'timestamp' column to datetime format.
- Computing the counting times and filling 'NaN' values in the counts for each detector.
- Normalizing the counts to counts per hour.
- Calculating the total counts across all detectors and discarding outliers.
- Finding similar stations based on cutoff rigidity, which is estimated using latitude and longitude values.
- Correcting the count by atmospheric variables and incoming neutron flux.
- Converting the smoothed neutron counts to Volumetric Water Content (VWC).

## [Rover Example: Hydroinnova Rover](../rover/Hydroinnova_rover_example/)

This example demonstrates how to process and analyze data from a rover using the Cosmic Ray Neutron Python (CRNPy) library. The steps include loading the data, converting geographical coordinates to Cartesian coordinates, normalizing neutron counts, correcting for atmospheric effects, and visualizing the results.

Key steps include:

- Data Loading: The data is loaded from a .csv file into a pandas DataFrame.

- Conversion to UTM Coordinates: The latitude and longitude are converted to UTM coordinates (x and y) for better plotting and analysis.

- Normalization of Neutron Counts: The neutron counts are normalized to counts per minute in case some observations covered a different timespan.

- Atmospheric Correction: The neutron counts are corrected for atmospheric variables and the incoming neutron flux.

- Smoothing and Conversion to VWC: The corrected counts are smoothed using a 2D smoothing function. The smoothed counts are then converted to volumetric water content (VWC) using the counts_to_vwc function.

- Interpolation and Visualization: The VWC is interpolated and a gridded map of the VWC is created using the interpolated x, y, and VWC values. The map is displayed using a color map to represent the VWC values. Finally, a contour map of the VWC is created and displayed.
