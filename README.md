# CRNPy
A Python toolbox for handling common tasks with cosmic-ray neutron probes. 

The toolbox was developed using common libraries for scientific computing like Pandas, NumPy, Matplotlib, and SciPy.

## Functionality
The toolbox seamlessly integrates with Pandas DataFrames for easy tabular data handling on top of the raw sensor data. The toolbox is desgined to read raw sensor data and append additional variables to the existing DataFrame. This way, researchers can export processed files while retaining all the raw data for better reproducibility and transparency.

- Helper functions for reading tabular data
- Fill timestamps with missing neutron counts
- Filter time series of neutron counts
- Corrections of neutron counts by atmospheric conditions and incoming neutron flux
- Conversion of neutron counts to volumetric water content
- Estimate sensing depth
- Estimate surface and sub-surface soil water storage
- Functions for calibration and validation of new detectors
- Out-of-the-box plotting functions to generate publication-quality figures


## References
These references represent peer-reviewed articles that used one or more functions included in the toolbox.

Dong, J., Ochsner, T.E., Zreda, M., Cosh, M.H. and Zou, C.B., 2014. Calibration and validation of the COSMOS rover for surface soil moisture measurement. Vadose Zone Journal, 13(4). doi.org/10.2136/vzj2013.08.0148

Patrignani, A., Ochsner, T., Montag, B. and Bellinger, S., 2021. A Novel Lithium Foil Cosmic-Ray Neutron Detector for Measuring Field-Scale Soil Moisture. Frontiers in Water, 3, p.67. doi.org/10.3389/frwa.2021.673185

Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020. Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9. doi.org/10.3389/frwa.2020.00009