![CRNPy Processing Workflow](img/workflow.png)
Overview of the proposed CRNPy processing workflow. Final user can choose to use the entire workflow, part of it, or build functions on top of it depending on their needs, dashed lines indicate optional steps.

-   Incoming neutron flux: The crnPy package includes a complete set of methods for correcting the raw observed neutron counts by natural variation in the incoming neutron flux, including a set of tools for searching and downloading data from a reference neutron monitoring station from the NMDB database (www.nmdb.eu) by proposing the most similar stations after analyzing the cut-off rigidity of the reference station and the estimated cut-off rigidity value for the studied location as a form of finding stations under a similar earth electromagnetic field (Klein et al., 2009; Shea & Smart, 2019; Smart & Shea, 2001).

| Incoming neutron flux correction factor|
|---------------------------------|
|$fi = \frac{I_{m}}{I_{0}}$|
|$ fi $: Incoming neutron flux correction factor|
|$ I_{m} $: Measured incoming neutron flux|
|$ I_{0} $: Reference incoming neutron flux at a given time.|

-   Atmospheric corrections: The package also provides functions for correcting raw neutron counts for atmospheric pressure, humidity, and temperature variations, (Andreasen et al., 2017; Rosolem et al., 2013; Zreda et al., 2012). Extra routines in this module account for other hydrogen environmental pools, lattice water, and total soil carbon that could impact cosmic-ray neutrons attenuation. A function to estimate the soil lattice water content (`estimate_lattice_water()`) based on texture and soil carbon was developed by analyzing collected soil samples across the state of Kansas combined with data published by Dong & Ochsner (2018)

| Pressure correction | Atmospheric water correction |
|---------------------|------------------------------|
|$fp = exp(\frac{P_{0} - P}{L})$ | $fw = 1 + 0.0054*(A - Aref)$ |
|$fp$: Atmospheric pressure correction factor | $fw$: Atmospheric water correction factor
|$P_{0}$: Reference atmospheric pressure (for e.g. long-term average) | $A$: Atmospheric water content
|$P$: Measured atmospheric pressure | $Aref$: Reference atmospheric water content
|$L$: Mass attenuation factor for high-energy neutrons in air | |

-   Biomass correction: The library provides a function for correcting neutron counts for the effects of above-ground biomass by combining an approach for estimating biomass water equivalent (BWE) from in-situ biomass samples (Wahbi et al., 2018) and the BWE correction factor described in Baatz et al. (2015).

| Biomass correction |
|--------------------|
|$fb = 1 - bwe*r2_N0$ |
|$fb$: Biomass correction factor |
|$bwe$: Biomass water equivalent |
|$r2_N0$: Ratio of the neutron counts reduction (counts kg-1) to the neutron calibration constant (N0) |

-   Road correction: Due to the known high sensibility closer to the detector, use cases like rover surveys could require the use of a correction factor accounting for the differences between the field soil water content and the road water content, crnPy implements the methodology proposed by Schrön et al. (2018).

| Road correction |
|-----------------|
|$fr = 1 + F1*F2*F3$ |
|$fr$: Road correction factor |
|$F1$: Road geometry term |
|$F2$: Road moisture term |
|$F3$: Road distance term |

See the original paper for more details, and the functions [crnpy.crnpy.road_correction][] documentation for the implementation details.

References:

Klein, K.-L., Steigies, C., & Nmdb Team. (2009). WWW.NMDB.EU: The real-time Neutron Monitor database. EGU General Assembly Conference Abstracts, 5633.

Shea, M., & Smart, D. (2019). Re-examination of the first five ground-level events. International Cosmic Ray Conference (ICRC2019), 36, 1149.

Smart, D., & Shea, M. (2001). Geomagnetic cutoff rigidity computer program: Theory, software description and example.

Andreasen, M., Jensen, K. H., Desilets, D., Franz, T. E., Zreda, M., Bogena, H. R., & Looms, M. C. (2017). Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone Journal, 16(8), 1–11.

Rosolem, R., Shuttleworth, W. J., Zreda, M., Franz, T. E., Zeng, X., & Kurc, S. A. (2013). The effect of atmospheric water vapor on neutron count in the cosmic-ray soil moisture observing system. Journal of Hydrometeorology, 14(5), 1659–1671.

Zreda, M., Desilets, D., Ferré, T., & Scott, R. L. (2008). Measuring soil moisture content non-invasively at intermediate spatial scale using cosmic-ray neutrons. Geophysical Research Letters, 35(21).

Dong, J., & Ochsner, T. E. (2018). Soil texture often exerts a stronger influence than precipitation on mesoscale soil moisture patterns. Water Resources Research, 54(3), 2199–2211.

Wahbi, A., Heng, L., Dercon, G., Wahbi, A., & Avery, W. (2018). In situ destructive sampling. Cosmic Ray Neutron Sensing: Estimation of Agricultural Crop Biomass Water Equivalent, 5–9.

Baatz, R., Bogena, H., Hendricks Franssen, H.-J., Huisman, J., Montzka, C., & Vereecken, H. (2015). An empirical vegetation correction for soil water content quantification using cosmic ray probes. Water Resources Research, 51(4), 2030–2046.

Schrön, M., Rosolem, R., Köhli, M., Piussi, L., Schröter, I., Iwema, J., Kögler, S., Oswald, S. E., Wollschläger, U., Samaniego, L., & others. (2018). Cosmic-ray neutron rover surveys of field soil moisture and the influence of roads. Water Resources Research, 54(9), 6441–6459.

Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., and Rosolem, R.: COSMOS: the COsmic-ray Soil Moisture Observing System, Hydrol. Earth Syst. Sci., 16, 4079–4099, https://doi.org/10.5194/hess-16-4079-2012, 2012.
