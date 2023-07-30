# crnpy/crnpy.py
"""
`crnpy` is a Python package for processing cosmic ray neutron data.

 Created by Joaquin Peraza and Andres Patrignani.
"""

# Import modules
import sys
import warnings
import numbers
import numpy as np
import pandas as pd
import requests
import io


from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from scipy.special import erfcinv
import crnpy.data as data

# Define python version
python_version = (3, 7)  # tuple of (major, minor) version requirement
python_version_str = str(python_version[0]) + "." + str(python_version[1])

# produce an error message if the python version is less than required
if sys.version_info < python_version:
    msg = "Module only runs on python version >= %s" % python_version_str
    raise Exception(msg)


def remove_incomplete_intervals(df, timestamp_col, integration_time, remove_first=False):
    """Function that removes rows with incomplete integration intervals.
    
    Args:
        df (pandas.DataFrame): Pandas Dataframe with data from stationary or roving CRNP devices.
        timestamp_col (str): Name of the column with timestamps in datetime format.
        integration_time (int): Duration of the neutron counting interval in seconds. Typical values are 60 seconds and 3600 seconds.
        remove_first (bool, optional): Remove first row. Default is False.
        
    Returns:
        (pandas.DataFrame): 
    """
    
    # Check format of timestamp column
    if df[timestamp_col].dtype != 'datetime64[ns]':
        raise TypeError('timestamp_col must be datetime64. Use `pd.to_datetime()` to fix this issue.')

    # Check if differences in timestamps are below or above the provided integration time
    idx_delta = df[timestamp_col].diff().dt.total_seconds() != integration_time
    
    if remove_first:
        idx_delta[0] = True
        
    # Select rows that meet the specified integration time
    df = df[~idx_delta]
    df.reset_index(drop=True, inplace=True)

    # Notify user about the number of rows that have been removed
    print(f"Removed a total of {sum(idx_delta)} rows.")
    
    return df


def fill_missing_timestamps(df, timestamp_col='timestamp', freq='H', round_timestamp=True, verbose=False):
    """Helper function to fill rows with missing timestamps in datetime record. Rows are filled with NaN values.

     Args:
         df (pandas.DataFrame): Pandas DataFrame.
         timestamp_col (str, optional): Column with the timestamp. Must be in datetime format. Default column name is 'timestamp'.
         freq (str, optional): Timestamp frequency. 'H' for hourly, 'M' for minute, or None. Can also use '3H' for a 3 hour frequency. Default is 'H'.
         round_timestamp (bool, optional): Whether to round timestamps to the nearest frequency. Default is True.
         verbose (bool, optional): Prints the missing timestamps added to the DatFrame.

     Returns:
         (pandas.DataFrame): DataFrame with filled missing timestamps.

     """

    # Check format of timestamp column
    if df[timestamp_col].dtype != 'datetime64[ns]':
        raise TypeError('timestamp_col must be datetime64. Use `pd.to_datetime()` to fix this issue.')

    # Round timestamps to nearest frequency. This steps must preced the filling of rows.
    if round_timestamp:
        df[timestamp_col] = df[timestamp_col].dt.round(freq)

    # Fill in rows with missing timestamps
    start_date = df[timestamp_col].iloc[0]
    end_date = df[timestamp_col].iloc[-1]
    date_range = pd.date_range(start_date, end_date, freq=freq)
    counter = 0
    for date in date_range:
        if date not in df[timestamp_col].values:
            if verbose:
                print('Adding missing date:',date)
            new_line = pd.DataFrame({timestamp_col:date}, index=[-1]) # By default fills columns with np.nan
            df = pd.concat([df,new_line])
            counter += 1

    df.sort_values(by=timestamp_col, inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # Notify user about the number of rows that have been removed
    print(f"Added a total of {counter} missing timestamps.")
    
    return df


def total_raw_counts(counts, nan_strategy=None, timestamp_col=None):
    """Compute the sum of uncorrected neutron counts for all detectors.

    Args:
        counts (pandas.DataFrame): Dataframe containing only the columns with neutron counts.
        nan_strategy (str): Strategy to use for NaN values. Options are 'interpolate', 'average', or None. Default is None.

    Returns:
        (pandas.DataFrame): Dataframe with the sum of uncorrected neutron counts for all detectors.
    """

    if counts.shape[0] > 1:
        counts = counts.apply(lambda x: x.fillna(counts.mean(axis=1)),axis=0)

    # Compute sum of counts
    total_raw_counts = counts.sum(axis=1)
    
    # Replace zeros with NaN
    total_raw_counts = total_raw_counts.replace(0, np.nan)
    
    return total_raw_counts


def is_outlier(x, method, window=11, min_val=None, max_val=None):
    """Function that tests whether values are outliers using a modified moving z-score based on the median absolute difference.

    Args:
        x (pandas.DataFrame): Variable containing only the columns with neutron counts.
        method (str): Outlier detection method. One of: range, iqr, moviqr, zscore, movzscore, modified_zscore, and scaled_mad
        window (int, optional): Window size for the moving central tendency. Default is 11.
        min_val (int or float): Minimum value for a reading to be considered valid. Default is None.
        max_val(int or float): Maximum value for a reading to be considered valid. Default is None.

    Returns:
        (pandas.DataFrame): Boolean indicating outliers.

    References:
        Iglewicz, B. and Hoaglin, D.C., 1993. How to detect and handle outliers (Vol. 16). Asq Press.
    """
        
    if not isinstance(x, pd.Series):
        raise TypeError('x must of type pandas.Series')

    # Separate this method to allow usage together with other methods below
    if isinstance(min_val, numbers.Number) and isinstance(max_val, numbers.Number):
        idx_range_outliers = (x < min_val) | (x > max_val)
    else:
        idx_range_outliers = np.full_like(x, False)

    # Apply other methods in addition to a range check
    if method == 'iqr':
        q1 = x.quantile(0.25)
        q3 = x.quantile(0.75)
        iqr = q3 - q1
        high_fence = q3 + (1.5 * iqr)
        low_fence = q1 - (1.5 * iqr)
        idx_outliers = (x<low_fence ) | (x>high_fence )

    elif method == 'moviqr':
        q1 = x.rolling(window, center=True).quantile(0.25)
        q3 = x.rolling(window, center=True).quantile(0.75)
        iqr = q3 - q1
        ub = q3 + (1.5 * iqr) # Upper boundary
        lb = q1 - (1.5 * iqr) # Lower boundary
        idx_outliers = (x < lb) | (x > ub)

    elif method == 'zscore':
        zscore = (x - x.mean())/x.std()
        idx_outliers = (zscore < -3) | (zscore > 3)

    elif method == 'movzscore':
        movmean = x.rolling(window=window, center=True).mean()
        movstd = x.rolling(window=window, center=True).std()
        movzscore = (x - movmean)/movstd
        idx_outliers = (movzscore < -3) | (movzscore > 3)

    elif method == 'modified_zscore':
        # Compute median absolute difference
        movmedian = x.rolling(window, center=True).median()
        abs_diff = np.abs(x - movmedian)
        mad = abs_diff.rolling(window, center=True).median()

        # Compute modified z-score
        modified_z_score = 0.6745 * abs_diff / mad
        idx_outliers = (modified_z_score < -3.5) | (modified_z_score > 3.5)

    elif method == 'scaled_mad':
        # Returns true for elements more than three scaled MAD from the median. 
        c = -1 / (np.sqrt(2)*erfcinv(3/2))
        median = np.nanmedian(x)
        mad = c*np.nanmedian(np.abs(x - median))
        idx_outliers = x > (median + 3*mad)

    else:
        raise TypeError('Outlier detection method not found.')

    return idx_outliers | idx_range_outliers


def correction_pressure(pressure, Pref, L):
    r"""Correction factor for atmospheric pressure.

    This function corrects neutron counts for atmospheric pressure using the method described in Andreasen et al. (2017).
    The correction is performed using the following equation:

    $$
    C_{corrected} = \frac{C_{raw}}{fp}
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fp: pressure correction factor

    $$
    fp = e^{\frac{P_{ref} - P}{L}}
    $$

    where:

    - P: atmospheric pressure
    - Pref: reference atmospheric pressure
    - L: Atmospheric attenuation coefficient.


    Args:
        atm_pressure (list or array): Atmospheric pressure readings. Long-term average pressure is recommended.
        Pref (float): Reference atmospheric pressure.
        L (float): Atmospheric attenuation coefficient.

    Returns:
        (list): fp pressure correction factor.

    References:
        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086
    """

    # Compute pressure correction factor
    fp = np.exp((Pref - pressure) / L) # Zreda et al. 2017 Eq 5.

    return fp


def correction_humidity(abs_humidity, Aref):
    r"""Correction factor for absolute humidity.

    This function corrects neutron counts for absolute humidity using the method described in Rosolem et al. (2013) and Anderson et al. (2017). The correction is performed using the following equation:

    $$
    C_{corrected} = C_{raw} \cdot f_w
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fw: absolute humidity correction factor

    $$
    f_w = 1 + 0.0054(A - A_{ref})
    $$

    where:

    - A: absolute humidity
    - Aref: reference absolute humidity

    Args:
        abs_humidity (list or array): Relative humidity readings.
        temp (list or array): Temperature readings (Celsius).
        Aref (float): Reference absolute humidity (g/m^3). The day of the instrument calibration is recommended.

    Returns:
        (list): fw correction factor.

    References:
        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086
    """
    A = abs_humidity
    fw = 1 + 0.0054*(A - Aref) # Zreda et al. 2017 Eq 6.
    return fw

def correction_incoming_flux(incoming_neutrons, incoming_Ref=None):
    r"""Correction factor for incoming neutron flux.

    This function corrects neutron counts for incoming neutron flux using the method described in Anderson et al. (2017). The correction is performed using the following equation:

    $$
    C_{corrected} = \frac{C_{raw}}{f_i}
    $$

    where:

    - Ccorrected: corrected neutron counts
    - Craw: raw neutron counts
    - fi: incoming neutron flux correction factor

    $$
    f_i = \frac{I_{ref}}{I}
    $$

    where:

    - I: incoming neutron flux
    - Iref: reference incoming neutron flux

    Args:
        incoming_neutrons (list or array): Incoming neutron flux readings.
        incoming_Ref (float): Reference incoming neutron flux. Baseline incoming neutron flux.

    Returns:
        (list): fi correction factor.

    References:
        M. Andreasen, K.H. Jensen, D. Desilets, T.E. Franz, M. Zreda, H.R. Bogena, and M.C. Looms. 2017. Status and perspectives on the cosmic-ray neutron method for soil moisture estimation and other environmental science applications. Vadose Zone J. 16(8). doi:10.2136/vzj2017.04.0086


    """
    if incoming_Ref is None and not isinstance(incoming_neutrons, type(None)):
        incoming_Ref = incoming_neutrons[0]
        warnings.warn('Reference incoming neutron flux not provided. Using first value of incoming neutron flux.')
    fi = incoming_neutrons / incoming_Ref
    fi.fillna(1.0, inplace=True)  # Use a value of 1 for days without data

    return fi


def get_incoming_neutron_flux(start_date, end_date, station, utc_offset=0, expand_window=0, verbose=False):
    """Function to retrieve neutron flux from the Neutron Monitor Database.

    Args:
        start_date (datetime): Start date of the time series.
        end_date (datetime): End date of the time series.
        station (str): Neutron Monitor station to retrieve data from.
        utc_offset (int): UTC offset in hours. Default is 0.
        expand_window (int): Number of hours to expand the time window to retrieve extra data. Default is 0.
        verbose (bool): Print information about the request. Default is False.

    Returns:
        (pandas.DataFrame): Neutron flux in counts per hour and timestamps.

    References:
        Documentation available:https://www.nmdb.eu/nest/help.php#howto
    """

    # Example: get_incoming_flux(station='IRKT',start_date='2020-04-10 11:00:00',end_date='2020-06-18 17:00:00')
    # Template url = 'http://nest.nmdb.eu/draw_graph.php?formchk=1&stations[]=KERG&output=ascii&tabchoice=revori&dtype=corr_for_efficiency&date_choice=bydate&start_year=2009&start_month=09&start_day=01&start_hour=00&start_min=00&end_year=2009&end_month=09&end_day=05&end_hour=23&end_min=59&yunits=0'


    # Expand the time window by 1 hour to ensure an extra observation is included in the request.
    start_date -= pd.Timedelta(hours=expand_window)
    end_date += pd.Timedelta(hours=expand_window)

    # Convert local time to UTC
    start_date = start_date - pd.Timedelta(hours=utc_offset)
    end_date = end_date - pd.Timedelta(hours=utc_offset)
    date_format = '%Y-%m-%d %H:%M:%S'
    root = 'http://www.nmdb.eu/nest/draw_graph.php?'
    url_par = [ 'formchk=1',
                'stations[]=' + station,
                'output=ascii',
                'tabchoice=revori',
                'dtype=corr_for_efficiency',
                'tresolution=' + str(60),
                'date_choice=bydate',
                'start_year=' + str(start_date.year),
                'start_month=' + str(start_date.month),
                'start_day=' + str(start_date.day),
                'start_hour=' + str(start_date.hour),
                'start_min=' + str(start_date.minute),
                'end_year=' + str(end_date.year),
                'end_month=' + str(end_date.month),
                'end_day=' + str(end_date.day),
                'end_hour=' + str(end_date.hour),
                'end_min=' + str(end_date.minute),
                'yunits=0']

    url = root + '&'.join(url_par)

    if verbose:
        print(f"Retrieving data from {url}")

    r = requests.get(url).content.decode('utf-8')

    # Subtract 1 hour to restore the last date included in the request.
    end_date -= pd.Timedelta('1H')
    start = r.find("RCORR_E\n") + 8
    end = r.find('\n</code></pre><br>Total') - 1
    s = r[start:end]
    s2 = ''.join([row.replace(';',',') for row in s])
    try:
        df_flux = pd.read_csv(io.StringIO(s2), names=['timestamp','counts'])
    except:
        if verbose:
            print(f"Error retrieving data from {url}")
        return None

    # Check if all values from selected detector are NaN. If yes, warn the user
    if df_flux['counts'].isna().all():
        warnings.warn('Data for selected neutron detectors appears to be unavailable for the selected period')

    # Convert timestamp to datetime and apply UTC offset
    df_flux['timestamp'] = pd.to_datetime(df_flux['timestamp'])
    df_flux['timestamp'] = df_flux['timestamp'] + pd.Timedelta(hours=utc_offset)

    # Print acknowledgement to inform users about restrictions and to acknowledge the NMDB database
    acknowledgement = """Data retrieved via NMDB are the property of the individual data providers. These data are free for non commercial
use to within the restriction imposed by the providers. If you use such data for your research or applications, please acknowledge
the origin by a sentence like 'We acknowledge the NMDB database (www.nmdb.eu) founded under the European Union's FP7 programme 
(contract no. 213007), and the PIs of individual neutron monitors at: IGY Jungfraujoch 
(Physikalisches Institut, University of Bern, Switzerland)"""

    return df_flux



def smooth_1d(values, window=5, order=3, method='moving_median'):
    """Use a Savitzky-Golay filter to smooth the signal of corrected neutron counts or another one-dimensional array (e.g. computed volumetric water content).

    Args:
        values (pd.DataFrame or pd.Serie): Dataframe containing the values to smooth.
        window (int): Window size for the Savitzky-Golay filter. Default is 5.
        method (str): Method to use for smoothing the data. Default is 'moving_median'.
            Options are 'moving_average', 'moving_median' and 'savitzky_golay'.
        order (int): Order of the Savitzky-Golay filter. Default is 3.

    Returns:
        (pd.DataFrame): DataFrame with smoothed values.

    References:
        Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020.
        Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.
        doi.org/10.3389/frwa.2020.00009
    """

    if method == 'moving_average':
        corrected_counts = values.rolling(window=window, center=True, min_periods=1).mean()
    elif method == 'moving_median':
        corrected_counts = values.rolling(window=window, center=True, min_periods=1).median()

    elif method == 'savitzky_golay':
        if values.isna().any():
            print('Dataframe contains NaN values. Please remove NaN values before smoothing the data.')

        if type(values) == pd.core.series.Series:
            filtered = savgol_filter(values,window,order)
            corrected_counts = pd.DataFrame(filtered,columns=['smoothed'], index=values.index)
        elif type(values) == pd.core.frame.DataFrame:
            for col in values.columns:
                values[col] = savgol_filter(values[col],window,order)
    else:
        raise ValueError('Invalid method. Please select a valid filtering method., options are: moving_average, moving_median, savitzky_golay')
    corrected_counts = corrected_counts.ffill(limit=window).bfill(limit=window).copy()
    return corrected_counts


def correction_bwe(counts, bwe, r2_N0=0.05):
    """Function to correct for biomass effects in neutron counts.
    following the approach described in Baatz et al., 2015.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of ephithermal neutron counts.
        bwe (float): Biomass water equivalent kg m-2.
        r2_N0 (float): Ratio of neutron counts with biomass to neutron counts without biomass. Default is 0.05.

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for biomass effects.

    References:
        Baatz, R., H. R. Bogena, H.-J. Hendricks Franssen, J. A. Huisman, C. Montzka, and H. Vereecken (2015),
        An empiricalvegetation correction for soil water content quantification using cosmic ray probes,
        Water Resour. Res., 51, 2030–2046, doi:10.1002/ 2014WR016443.
    """

    return counts/(1 - bwe*r2_N0)

def biomass_to_bwe(biomass_dry, biomass_fresh, fWE=0.494):
    """Function to convert biomass to biomass water equivalent.

    Args:
        biomass_dry (array or pd.Series or pd.DataFrame): Above ground dry biomass in kg m-2.
        biomass_fresh (array or pd.Series or pd.DataFrame): Above ground fresh biomass in kg m-2.
        fWE (float): Stoichiometric ratio of H2O to organic carbon molecules in the plant (assuming this is mostly cellulose)
            Default is 0.494 (Wahbi & Avery, 2018).

    Returns:
        (array or pd.Series or pd.DataFrame): Biomass water equivalent in kg m-2.

    References:
        Wahbi, A., Avery, W. (2018). In Situ Destructive Sampling. In:
        Cosmic Ray Neutron Sensing: Estimation of Agricultural Crop Biomass Water Equivalent.
        Springer, Cham. https://doi.org/10.1007/978-3-319-69539-6_2
    """
    return (biomass_fresh - biomass_dry) + fWE * biomass_dry


def correction_road(counts, theta_N, road_width, road_distance=0.0, theta_road=0.12, p0=0.42, p1=0.5, p2=1.06, p3=4, p4=0.16, p6=0.94, p7=1.10, p8=2.70, p9=0.01):
    """Function to correct for road effects in neutron counts.
    following the approach described in Schrön et al., 2018.

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of ephithermal neutron counts.
        theta_N (float): Volumetric water content of the soil estimated from the uncorrected neutron counts.
        road_width (float): Width of the road in m.
        road_distance (float): Distance of the road from the sensor in m. Default is 0.0.
        theta_road (float): Volumetric water content of the road. Default is 0.12.
        p0-p9 (float): Parameters of the correction function. Default values are from Schrön et al., 2018.

    Returns:
        (array or pd.Series or pd.DataFrame): Array of corrected neutron counts for road effects.

    References:
        Schrön,M.,Rosolem,R.,Köhli,M., Piussi,L.,Schröter,I.,Iwema,J.,etal. (2018).Cosmic-ray neutron rover surveys
        of field soil moisture and the influence of roads.WaterResources Research,54,6441–6459.
        https://doi. org/10.1029/2017WR021719
    """
    F1 = p0 * (1-np.exp(-p1*road_width))
    F2 = -p2 - p3 * theta_road - ((p4 + theta_road) / (theta_N))
    F3 = p6 * np.exp(-p7 * (road_width ** -p8) * road_distance ** 4) + (1 - p6) * np.exp(-p9 * road_distance)

    C_roads = 1 + F1 * F2 * F3

    corrected_counts = counts / C_roads

    return corrected_counts

def counts_to_vwc(counts, N0, Wlat, Wsoc ,bulk_density, a0=0.0808,a1=0.372,a2=0.115):
    r"""Function to convert corrected and filtered neutron counts into volumetric water content.

    This method implements soil moisture estimation using the non-linear relationship between neutron count and soil volumetric water content following the approach described in Desilets et al., 2010.

    $\theta(N) =\frac{a_0}{(\frac{N}{N_0}) - a_1} - a_2 $

    Args:
        counts (array or pd.Series or pd.DataFrame): Array of corrected and filtered neutron counts.
        N0 (float): Device-specific neutron calibration constant.
        Wlat (float): Lattice water content.
        Wsoc (float): Soil organic carbon content.
        bulk_density (float): Soil bulk density.
        a0 (float): Parameter given in Zreda et al., 2012. Default is 0.0808.
        a1 (float): Parameter given in Zreda et al., 2012. Default is 0.372.
        a2 (float): Parameter given in Zreda et al., 2012. Default is 0.115.

    Returns:
        (array or pd.Series or pd.DataFrame): Volumetric water content in m3 m-3.

    References:
        Desilets, D., M. Zreda, and T.P.A. Ferré. 2010. Nature’s neutron probe:
        Land surface hydrology at an elusive scale with cosmic rays. Water Resour. Res. 46:W11505.
        doi.org/10.1029/2009WR008726
    """

    # Convert neutron counts into vwc
    vwc = (a0 / (counts/N0-a1) - a2 - Wlat - Wsoc) * bulk_density
    return vwc



def sensing_depth(vwc, pressure, p_ref, bulk_density, Wlat, dist=None, method='Schron_2017'):
    # Convert docstring to google format
    """Function that computes the estimated sensing depth of the cosmic-ray neutron probe.
    The function offers several methods to compute the depth at which 86 % of the neutrons
    probe the soil profile.

    Args:
        vwc (array or pd.Series or pd.DataFrame): Estimated volumetric water content for each timestamp.
        pressure (array or pd.Series or pd.DataFrame): Atmospheric pressure in hPa for each timestamp.
        p_ref (float): Reference pressure in hPa.
        bulk_density (float): Soil bulk density.
        Wlat (float): Lattice water content.
        method (str): Method to compute the sensing depth. Options are 'Schron_2017' or 'Franz_2012'.
        dist (list or array): List of radial distances at which to estimate the sensing depth. Only used for the 'Schron_2017' method.

    Returns:
        (array or pd.Series or pd.DataFrame): Estimated sensing depth in m.

    References:
        Franz, T.E., Zreda, M., Ferre, T.P.A., Rosolem, R., Zweck, C., Stillman, S., Zeng, X. and Shuttleworth, W.J., 2012.
        Measurement depth of the cosmic ray soil moisture probe affected by hydrogen from various sources.
        Water Resources Research, 48(8). doi.org/10.1029/2012WR011871

        Schrön, M., Köhli, M., Scheiffele, L., Iwema, J., Bogena, H. R., Lv, L., et al. (2017).
        Improving calibration and validation of cosmic-ray neutron sensors in the light of spatial sensitivity.
        Hydrol. Earth Syst. Sci. 21, 5009–5030. doi.org/10.5194/hess-21-5009-2017
    """

    # Determine sensing depth (D86)
    if method == 'Schron_2017':
        # See Appendix A of Schrön et al. (2017)
        Fp = 0.4922 / (0.86 - np.exp(-1 * pressure / p_ref));
        Fveg = 0
        results = []
        for d in dist:
            # Compute r_star
            r_start = d/Fp

            # Compute soil depth that accounts for 86% of the neutron flux
            D86 = 1/ bulk_density * (8.321+0.14249*(0.96655 + np.exp(-0.01*r_start))*(20+(Wlat+vwc)) / (0.0429+(Wlat+vwc)))
            results.append(D86)

    elif method == 'Franz_2012':
        results = 5.8/(bulk_density*Wlat+vwc+0.0829)
    else:
        raise ValueError('Method not recognized. Please select either "Schron_2017" or "Franz_2012".')

    return results

def abs_humidity(relative_humidity, temp):
    """
    Compute the actual vapor pressure (e) in g m^-3 using RH (%) and current temperature (c) observations.

    Args:
        relative_humidity (float): relative humidity (%)
        temp (float): temperature (Celsius)

    Returns:
        float: actual vapor pressure (g m^-3)
    """

    ### Atmospheric water vapor factor
    # Saturation vapor pressure
    e_sat = 0.611 * np.exp(17.502 * temp / (
                temp + 240.97)) * 1000  # in Pascals Eq. 3.8 p.41 Environmental Biophysics (Campbell and Norman)

    # Vapor pressure Pascals
    Pw = e_sat * relative_humidity / 100

    # Absolute humidity (g/m^3)
    C = 2.16679  # g K/J;
    abs_h = C * Pw / (temp + 273.15)
    return abs_h


def nrad_weight(h,theta,distances,depth,rhob=1.4):
    """Function to compute distance weights corresponding to each soil sample.

    Args:
        h (float): Air Humidity  from 0.1  to 50 in g/m^3. When h=0, the function will skip the distance weighting.
        theta (array or pd.Series or pd.DataFrame): Soil Moisture for each sample (0.02 - 0.50 m^3/m^3)
        distances (array or pd.Series or pd.DataFrame): Distances from the location of each sample to the origin (0.5 - 600 m)
        depth (array or pd.Series or pd.DataFrame): Depths for each sample (m)
        rhob (float): Bulk density in g/cm^3

    Returns:
        (array or pd.Series or pd.DataFrame): Distance weights for each sample.

    References:
        Köhli, M., Schrön, M., Zreda, M., Schmidt, U., Dietrich, P., and Zacharias, S. (2015).
        Footprint characteristics revised for field-scale soil moisture monitoring with cosmic-ray
        neutrons. Water Resour. Res. 51, 5772–5790. doi:10.1002/2015WR017169
    """

    # Table A1. Parameters for Fi and D86
    p10 = 8735;       p11 = 17.1758; p12 = 11720;      p13 = 0.00978;   p14 = 7045;      p15 = 0.003632;
    p20 = 2.7925e-2;  p21 = 5.0399;  p22 = 2.8544e-2;  p23 = 0.002455;  p24 = 6.851e-5;  p25 = 9.2926;
    p30 = 247970;     p31 = 17.63;   p32 = 374655;     p33 = 0.00191;   p34 = 195725;
    p40 = 5.4818e-2;  p41 = 15.921;  p42 = 0.6373;     p43 = 5.99e-2;   p44 = 5.425e-4;
    p50 = 1383702;    p51 = 4.156;   p52 = 5325;       p53 = 0.00238;   p54 = 0.0156;    p55 = 0.130;     p56 = 1521;
    p60 = 6.031e-5;   p61 = 98.5;    p62 = 1.0466e-3;
    p70 = 11747;      p71 = 41.66;   p72 = 4521;       p73 = 0.01998;   p74 = 0.00604;   p75 = 2534;      p76 = 0.00475;
    p80 = 1.543e-2;   p81 = 10.06;   p82 = 1.807e-2;   p83 = 0.0011;    p84 = 8.81e-5;   p85 = 0.0405;    p86 = 20.24;
    p90 = 8.321;      p91 = 0.14249; p92 = 0.96655;    p93 = 26.42;     p94 = 0.0567;


    # Numerical determination of the penetration depth (86%) (Eq. 8)
    D86 = 1/rhob*(p90+p91*(p92+np.exp(-1*distances/100))*(p93+theta)/(p94+theta))

    # Depth weights (Eq. 7)
    Wd = np.exp(-2*depth/D86)

    if h == 0:
        W = 1 # skip distance weighting

    elif (h >= 0.1) and (h<= 50):
        # Functions for Fi (Appendix A in Köhli et al., 2015)
        F1 = p10*(1+p13*h)*np.exp(-p11*theta)+p12*(1+p15*h)-p14*theta
        F2 = ((-p20+p24*h)*np.exp(-p21*theta/(1+p25*theta))+p22)*(1+h*p23)
        F3 = (p30*(1+p33*h)*np.exp(-p31*theta)+p32-p34*theta)
        F4 = p40*np.exp(-p41*theta)+p42-p43*theta+p44*h
        F5 = p50*(0.02-1/p55/(h-p55+p56*theta))*(p54-theta)*np.exp(-p51*(theta-p54))+p52*(0.7-h*theta*p53)
        F6 = p60*(h+p61)+p62*theta
        F7 = (p70*(1-p76*h)*np.exp(-p71*theta*(1-h*p74))+p72-p75*theta)*(2+h*p73)
        F8 = ((-p80+p84*h)*np.exp(-p81*theta/(1+p85*h+p86*theta))+p82)*(2+h*p83)

        # Distance weights (Eq. 3)
        W = np.ones_like(distances)*np.nan
        for i in range(len(distances)):
            if (distances[i]<=50) and (distances[i]>0.5):
                W[i]=F1[i]*(np.exp(-F2[i]*distances[i]))+F3[i]*np.exp(-F4[i]*distances[i])

            elif (distances[i]>50) and (distances[i]<600):
                W[i]=F5[i]*(np.exp(-F6[i]*distances[i]))+F7[i]*np.exp(-F8[i]*distances[i])

            else:
                raise ValueError('Input distances are not valid.')

    else:
        raise ValueError('Air humidity values are out of range.')


    # Combined and normalized weights
    weights = Wd*W/np.nansum(Wd*W)

    return weights



def exp_filter(sm,T=1):
    """Exponential filter to estimate soil moisture in the rootzone from surface observtions.

    Args:
        sm (list or array): Soil moisture in mm of water for the top layer of the soil profile.
        T (float): Characteristic time length in the same units as the measurement interval.

    Returns:
        sm_subsurface (list or array): Subsurface soil moisture in the same units as the input.

    References:
        Albergel, C., Rüdiger, C., Pellarin, T., Calvet, J.C., Fritz, N., Froissard, F., Suquia, D., Petitpa, A., Piguet, B. and Martin, E., 2008.
        From near-surface to root-zone soil moisture using an exponential filter: an assessment of the method based on in-situ observations and model
        simulations. Hydrology and Earth System Sciences, 12(6), pp.1323-1337.

        Franz, T.E., Wahbi, A., Zhang, J., Vreugdenhil, M., Heng, L., Dercon, G., Strauss, P., Brocca, L. and Wagner, W., 2020.
        Practical data products from cosmic-ray neutron sensing for hydrological applications. Frontiers in Water, 2, p.9.

        Rossini, P. and Patrignani, A., 2021. Predicting rootzone soil moisture from surface observations in cropland using an exponential filter.
        Soil Science Society of America Journal.
    """


    # Parameters
    t_delta = 1
    sm_min = np.min(sm)
    sm_max = np.max(sm)
    ms = (sm - sm_min) / (sm_max - sm_min)

    # Pre-allocate soil water index array and recursive constant K
    SWI = np.ones_like(ms)*np.nan
    K = np.ones_like(ms)*np.nan

    # Initial conditions
    SWI[0] = ms[0]
    K[0] = 1

    # Values from 2 to N
    for n in range(1,len(SWI)):
        if ~np.isnan(ms[n]) & ~np.isnan(ms[n-1]):
            K[n] = K[n-1] / (K[n-1] + np.exp(-t_delta/T))
            SWI[n] = SWI[n-1] + K[n]*(ms[n] - SWI[n-1])
        else:
            continue

    # Rootzone storage
    sm_subsurface = SWI * (sm_max - sm_min) + sm_min

    return sm_subsurface


def cutoff_rigidity(lat,lon):
    """Function to estimate the approximate cutoff rigidity for any point on Earth according to the
    tabulated data of Smart and Shea, 2019. The returned value can be used to select the appropriate
    neutron monitor station to estimate the cosmic-ray neutron intensity at the location of interest.

    Args:
        lat (float): Geographic latitude in decimal degrees. Value in range -90 to 90
        lon (float): Geographic longitude in decimal degrees. Values in range from 0 to 360.
            Typical negative longitudes in the west hemisphere will fall in the range 180 to 360.

    Returns:
        (float): Cutoff rigidity in GV. Error is about +/- 0.3 GV

    Examples:
        Estimate the cutoff rigidity for Newark, NJ, US

        >>> zq = cutoff_rigidity(39.68, -75.75)
        >>> print(zq)
        2.52 GV (Value from NMD is 2.40 GV)

    References:
        Smart, D. & Shea, Matthew. (2001). Geomagnetic Cutoff Rigidity Computer Program:
        Theory, Software Description and Example. NASA STI/Recon Technical Report N.

        Shea, M. A., & Smart, D. F. (2019, July). Re-examination of the First Five Ground-Level Events.
        In International Cosmic Ray Conference (ICRC2019) (Vol. 36, p. 1149).
    """
    xq = lon
    yq = lat

    if xq < 0:
        xq = xq*-1 + 180
    Z = np.array(data.cutoff_rigidity)
    x = np.linspace(0, 360, Z.shape[1])
    y = np.linspace(90, -90, Z.shape[0])
    X, Y = np.meshgrid(x, y)
    points = np.array( (X.flatten(), Y.flatten()) ).T
    values = Z.flatten()
    zq = griddata(points, values, (xq,yq))

    return np.round(zq,2)


def find_neutron_monitor(Rc, start_date=None, end_date=None, verbose=False):
    """Search for potential reference neutron monitoring stations based on cutoff rigidity.

    Args:
        Rc (float): Cutoff rigidity in GV. Values in range 1.0 to 3.0 GV.
        start_date (datetime): Start date for the period of interest.
        end_date (datetime): End date for the period of interest.

    Returns:
        (list): List of top five stations with closes cutoff rigidity.
            User needs to select station according to site altitude.

    Examples:
        >>> from crnpy import crnpy
        >>> Rc = 2.40 # 2.40 Newark, NJ, US
        >>> crnpy.find_neutron_monitor(Rc)
        Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations

        Your cutoff rigidity is 2.4 GV
                STID                          NAME     R  Altitude_m
        40   NEWK                        Newark  2.40          50
        33   MOSC                        Moscow  2.43         200
        27   KIEL                          Kiel  2.36          54
        28  KIEL2                        KielRT  2.36          54
        31   MCRL  Mobile Cosmic Ray Laboratory  2.46         200
        32   MGDN                       Magadan  2.10         220
        42   NVBK                   Novosibirsk  2.91         163
        26   KGSN                      Kingston  1.88          65
        9    CLMX                        Climax  3.00        3400
        57   YKTK                       Yakutsk  1.65         105

    References:
        https://www.nmdb.eu/nest/help.php#helpstations
    """

    # Load file with list of neutron monitoring stations
    stations = pd.DataFrame(data.neutron_detectors, columns=["STID","NAME","R","Altitude_m"])

    # Sort stations by closest cutoff rigidity
    idx_R = (stations['R'] - Rc).abs().argsort()

    if start_date is not None and end_date is not None:
        stations["Period available"] = False
        for i in range(10):
            station = stations.iloc[idx_R[i]]["STID"]
            try:
                if get_incoming_neutron_flux(start_date, end_date, station, verbose=verbose) is not None:
                    stations.iloc[idx_R[i],-1] = True
            except:
                pass
        if sum(stations["Period available"] == True) == 0:
            print("No stations available for the selected period!")
        else:
            stations = stations[stations["Period available"] == True]
            idx_R = (stations['R'] - Rc).abs().argsort()
            result = stations.iloc[idx_R.iloc[:10]]
    else:
        result = stations.reindex(idx_R).head(10).rename_axis(None)

    # Print results
    print('')
    print("""Select a station with an altitude similar to that of your location. For more information go to: 'https://www.nmdb.eu/nest/help.php#helpstations""")
    print('')
    print(f"Your cutoff rigidity is {Rc} GV")
    print(result)
    return result


def interpolate_incoming_flux(nmdb_timestamps, nmdb_counts, crnp_timestamps):
    """Function to interpolate incoming neutron flux to match the timestamps of the observations.

    Args:
        nmdb_timestamps (pd.Series): Series of timestamps in datetime format from the NMDB.
        nmdb_counts (pd.Series): Series of neutron counts from the NMDB
        crnp_timestamps (pd.Series): Series of timestamps in datetime format from the station or device.

    Returns:
        (pd.Series): Series containing interpolated incoming neutron flux. Length of Series is the same as crnp_timestamps
    """
    incoming_flux = np.array([])
    for k,timestamp in enumerate(crnp_timestamps):
        if timestamp in nmdb_timestamps.values:
            idx = timestamp == nmdb_timestamps
            incoming_flux = np.append(incoming_flux, nmdb_counts.loc[idx])
        else:
            incoming_flux = np.append(incoming_flux, np.nan)

    # Interpolate nan values
    incoming_flux = pd.Series(incoming_flux).interpolate(method='nearest', limit_direction='both')

    # Return only the values for the selected timestamps
    return incoming_flux


def lattice_water(clay_content, total_carbon=None):
    r"""Estimate the amount of water in the lattice of clay minerals.

    ![img1](img/lattice_water_simple.png) | ![img2](img/lattice_water_multiple.png)
    :-------------------------:|:-------------------------:
    $\omega_{lat} = 0.097 * clay(\%)$ | $\omega_{lat} = -0.028 + 0.077 * clay(\%) + 0.459 * carbon(\%)$
    Linear regression [lattice water (%) as a function of clay (%)] done with data from Kansas Sate University - Soil Water Processes Lab. |  Multiple linear regression [lattice water (%) as a function of clay (%) and soil carbon (%)] done with data from Soil Water Processes Lab.

    Args:
        clay_content (float): Clay content in the soil in percent.
        total_carbon (float, optional): Total carbon content in the soil in percent.
            If None, the amount of water is estimated based on clay content only.

    Returns:
        (float): Amount of water in the lattice of clay minerals in percent
    """
    if total_carbon is None:
        lattice_water = 0.097 * clay_content
    else:
        lattice_water = -0.028 + 0.077 * clay_content + 0.459 * total_carbon
    return lattice_water


def latlon_to_utm(lat, lon, utm_zone_number, missing_values=None):
    """Convert geographic coordinates (lat, lon) to projected coordinates (utm) using the Military Grid Reference System.

    Function only applies to non-polar coordinates.
    If further functionality is required, consider using the utm module. See references for more information.

    ![UTM zones](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b7/Universal_Transverse_Mercator_zones.svg/1920px-Universal_Transverse_Mercator_zones.svg.png)
    UTM zones on an equirectangular world map with irregular zones in red and New York City's zone highlighted. See [UTM zones](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#UTM_zones) for a full description.


    Args:
        lat (float, array): Latitude in decimal degrees.
        lon (float, array): Longitude in decimal degrees.
        utm_zone_number (int): Universal Transverse Mercator (UTM) zone.

    Returns:
        (float, float): Tuple of easting and northing coordinates in meters. First element is easting, second is northing.

    References:
         Code adapted from utm module created by Tobias Bieniek (Github username: Turbo87)
         [https://github.com/Turbo87/utm](https://github.com/Turbo87/utm)

         [https://www.maptools.com/tutorials/grid_zone_details#](https://www.maptools.com/tutorials/grid_zone_details#)
    """


    # Define constants
    R = 6_378_137  # Earth's radius at the Equator in meters

    # Convert input data to Numpy arrays
    if (type(lat) is not np.ndarray) or (type(lon) is not np.ndarray):
        try:
            lat = np.array(lat)
            lon = np.array(lon)
        except:
            raise "Input values cannot be converted to Numpy arrays."

    # Check latitude range
    if np.any(lat < -80) | np.any(lat > 84):
        raise "One or more latitude values exceed the range -80 to 84"

    # Check longitude range
    if np.any(lon < -180) | np.any(lon > 180):
        raise "One or more longitude values exceed the range -180 to 180"

    # Constants
    K0 = 0.9996
    E = 0.00669438
    E_P2 = E / (1 - E)

    M1 = (1 - E / 4 - 3 * E ** 2 / 64 - 5 * E ** 3 / 256)
    M2 = (3 * E / 8 + 3 * E ** 2 / 32 + 45 * E ** 3 / 1024)
    M3 = (15 * E ** 2 / 256 + 45 * E ** 3 / 1024)
    M4 = (35 * E ** 3 / 3072)

    # Trigonometric operations
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    lat_sin = np.sin(lat_rad)
    lat_cos = np.cos(lat_rad)
    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    # Find central meridian.
    central_lon = (utm_zone_number * 6 - 180) - 3  # Zones are every 6 degrees.
    central_lon_rad = np.radians(central_lon)

    n = R / np.sqrt(1 - E * lat_sin ** 2)
    c = E_P2 * lat_cos ** 2

    with np.errstate(divide='ignore', invalid='ignore'):
        a = lat_cos * (np.remainder(((lon_rad - central_lon_rad) + np.pi), (2 * np.pi)) - np.pi)
    m = R * (M1 * lat_rad - M2 * np.sin(2 * lat_rad) + M3 * np.sin(4 * lat_rad) - M4 * np.sin(6 * lat_rad))

    easting = K0 * n * (a + a ** 3 / 6 * (1 - lat_tan2 + c) + a ** 5 / 120 * (
                5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500_000
    northing = K0 * (m + n * lat_tan * (
                a ** 2 / 2 + a ** 4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c ** 2) + a ** 6 / 720 * (
                    61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))

    if np.any(lat < 0):
        northing += 10_000_000

    return easting, northing

def euclidean_distance(px, py, x, y):
    """Function that computes the Euclidean distance between one point
    in space and one or more points.

    Args:
        px (float): x projected coordinate of the point.
        py (float): y projected coordinate of the point.
        x (list, ndarray, pandas.series): vector of x projected coordinates.
        y (list, ndarray, pandas.series): vector of y projected coordinates.

    Returns:
        (ndarray): Numpy array of distances from the point (px,py) to all the points in x and y vectors.
    """
    d = np.sqrt((px - x) ** 2 + (py - y) ** 2)
    return d


def spatial_average(x, y, z, buffer=100, min_neighbours=3, method='mean', rnd=False):
    """Moving buffer filter to smooth georeferenced two-dimensional data.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be smoothed.
        buffer (float): Radial buffer distance in meters.
        min_neighbours (int): Minimum number of neighbours to consider for the smoothing.
        method (str): One of 'mean' or 'median'.
        rnd (bool): Boolean to round the final result. Useful in case of z representing neutron counts.

    Returns:
        (array): Smoothed version of z with the same dimension as z.
    """

    # Convert input data to Numpy arrays
    if (type(x) is not np.ndarray) or (type(y) is not np.ndarray):
        try:
            x = np.array(x)
            y = np.array(y)
        except:
            raise "Input values cannot be converted to Numpy arrays."

    if len(x) != len(y):
        raise f"The number of x and y must be equal. Input x has {len(x)} values and y has {len(y)} values."

    # Compute distances
    N = len(x)
    z_smooth = np.array([])
    for k in range(N):
        px = x[k]
        py = y[k]

        distances = euclidean_distance(px, py, x, y)
        idx_within_buffer = distances <= buffer


        if np.isnan(z[k]):
            z_new_val = np.nan
        elif len(distances[idx_within_buffer]) > min_neighbours:
            if method == 'mean':
                z_new_val = np.nanmean(z[idx_within_buffer])
            elif method == 'median':
                z_new_val = np.nanmedian(z[idx_within_buffer])
            else:
                raise f"Method {method} does not exist. Provide either 'mean' or 'median'."
        else:
            z_new_val = z[k] # If there are not enough neighbours, keep the original value

        # Append smoothed value to array
        z_smooth = np.append(z_smooth, z_new_val)

    if rnd:
        z_smooth = np.round(z_smooth, 0)

    return z_smooth


def idw(x, y, z, X_pred, Y_pred, neighborhood=1000, p=1):
    """Function to interpolate data using inverse distance weight.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be interpolated.
        X_pred (list or array): UTM x coordinates where z values need to be predicted.
        Y_pred (list or array): UTM y coordinates where z values need to be predicted.
        neighborhood (float): Only points within this radius in meters are considered for the interpolation.
        p (int): Exponent of the inverse distance weight formula. Typically, p=1 or p=2.

    Returns:
        (array): Interpolated values.

    References:
        [https://en.wikipedia.org/wiki/Inverse_distance_weighting](https://en.wikipedia.org/wiki/Inverse_distance_weighting)


    """

    # Flatten arrays to handle 1D and 2D arrays with the same code
    s = X_pred.shape  # Save shape
    X_pred = X_pred.flatten()
    Y_pred = Y_pred.flatten()

    # Pre-allocate output array
    Z_pred = np.full_like(X_pred, np.nan)

    for n in range(X_pred.size):
        # Distance between current and observed points
        d = euclidean_distance(X_pred[n], Y_pred[n], x, y)

        # Select points within neighborhood only for interpolation
        idx_neighbors = d < neighborhood

        # Compute interpolated value at point of interest
        Z_pred[n] = np.sum(z[idx_neighbors] / d[idx_neighbors] ** p) / np.sum(1 / d[idx_neighbors] ** p)

    return np.reshape(Z_pred, s)


def interpolate_2d(x, y, z, dx=100, dy=100, method='cubic', neighborhood=1000):
    """Function for interpolating irregular spatial data into a regular grid.

    Args:
        x (list or array): UTM x coordinates in meters.
        y (list or array): UTM y coordinates in meters.
        z (list or array): Values to be interpolated.
        dx (float): Pixel width in meters.
        dy (float): Pixel height in meters.
        method (str): Interpolation method. One of 'cubic', 'linear', 'nearest', or 'idw'.
        neighborhood (float): Only points within this radius in meters are considered for the interpolation.

    Returns:
        x_pred (array): 2D array with x coordinates.
        y_pred (array): 2D array with y coordinates.
        z_pred (array): 2D array with interpolated values.

    References:
        [https://soilwater.github.io/pynotes-agriscience/notebooks/interpolation.html](https://soilwater.github.io/pynotes-agriscience/notebooks/interpolation.html)
    """

    # Drop NaN values in x y and z
    idx_nan = np.isnan(x) | np.isnan(y) | np.isnan(z)
    x = x[~idx_nan]
    y = y[~idx_nan]
    z = z[~idx_nan]

    if idx_nan.any():
        print(f"WARNING: {np.isnan(x).sum()}, {np.isnan(y).sum()}, and {np.isnan(z).sum()} NaN values were dropped from x, y, and z.")

    # Create 2D grid for interpolation
    Nx = round((np.max(x) - np.min(x)) / dx) + 1
    Ny = round((np.max(y) - np.min(y)) / dy) + 1
    X_vec = np.linspace(np.min(x), np.max(x), Nx)
    Y_vec = np.linspace(np.min(y), np.max(y), Ny)
    X_pred, Y_pred = np.meshgrid(X_vec, Y_vec)

    if method in ['linear', 'nearest', 'cubic']:
        points = list(zip(x, y))
        Z_pred = griddata(points, z, (X_pred, Y_pred), method=method)

    elif method == 'idw':
        Z_pred = idw(x, y, z, X_pred, Y_pred, neighborhood)

    else:
        raise f"Method {method} does not exist. Provide either 'cubic', 'linear', 'nearest', or 'idw'."

    return X_pred, Y_pred, Z_pred


def rover_centered_coordinates(x, y):
    """Function to estimate the intermediate locations between two points, assuming the measurements were taken at a constant speed.

    Args:
        x (array): x coordinates.
        y (array): y coordinates.

    Returns:
        x_est (array): Estimated x coordinates.
        y_est (array): Estimated y coordinates.
    """

    # Make it datatype agnostic
    if(isinstance(x, pd.Series)):
        x = x.values
    if(isinstance(y, pd.Series)):
        y = y.values

    # Do the average of the two points
    x_est = (x[1:] + x[:-1]) / 2
    y_est = (y[1:] + y[:-1]) / 2

    # Add the first point to match the length of the original array
    x_est = np.insert(x_est, 0, x[0])
    y_est = np.insert(y_est, 0, y[0])


    return x_est, y_est


def uncertainty_counts(raw_counts, metric="std", fp=1, fw=1, fi=1):
    """Function to estimate the uncertainty of raw counts.

    Measurements of proportional neutron detector systems are governed by counting statistics that follow a Poissonian probability distribution (Zreda et al., 2012).
    The expected uncertainty in the neutron count rate $N$ is defined by the standard deviation $ \sqrt{N} $ (Jakobi et al., 2020).
    The CV% can be expressed as $ N^{-1/2} $

    Args:
        raw_counts (array): Raw neutron counts.
        metric (str): Either 'std' or 'cv' for standard deviation or coefficient of variation.
        fp (float): Pressure correction factor.
        fw (float): Humidity correction factor.
        fi (float): Incoming neutron flux correction factor.

    Returns:
        uncertainty (float): Uncertainty of raw counts.

    References:
        Jakobi J, Huisman JA, Schrön M, Fiedler J, Brogi C, Vereecken H and Bogena HR (2020) Error Estimation for Soil Moisture Measurements With
        Cosmic Ray Neutron Sensing and Implications for Rover Surveys. Front. Water 2:10. doi: 10.3389/frwa.2020.00010

        Zreda, M., Shuttleworth, W. J., Zeng, X., Zweck, C., Desilets, D., Franz, T., and Rosolem, R.: COSMOS: the COsmic-ray Soil Moisture Observing System,
        Hydrol. Earth Syst. Sci., 16, 4079–4099, https://doi.org/10.5194/hess-16-4079-2012, 2012.

    """

    s = fw / (fp * fi)
    if metric == "std":
        uncertainty = np.sqrt(raw_counts) * s
    elif metric == "cv":
        uncertainty = 1 /  np.sqrt(raw_counts) * s
    else:
        raise f"Metric {metric} does not exist. Provide either 'std' or 'cv' for standard deviation or coefficient of variation."
    return uncertainty


def uncertainty_vwc(raw_counts, N0, bulk_density, fp=1, fw=1, fi=1, a0=0.0808,a1=0.372,a2=0.115):
    r"""Function to estimate the uncertainty propagated to volumetric water content.

    The uncertainty of the volumetric water content is estimated by propagating the uncertainty of the raw counts.
    Following Eq. 10 in Jakobi et al. (2020), the uncertainty of the volumetric water content can be expressed as:
    $$
    \sigma_{\theta_g}(N) = \sigma_N \frac{a_0 N_0}{(N_{cor} - a_1 N_0)^4} \sqrt{(N_{cor} - a_1 N_0)^4 + 8 \sigma_N^2 (N_{cor} - a_1 N_0)^2 + 15 \sigma_N^4}
    $$

    Args:
        raw_counts (array): Raw neutron counts.
        N0 (float): Calibration parameter N0.
        bulk_density (float): Bulk density in g cm-3.
        fp (float): Pressure correction factor.
        fw (float): Humidity correction factor.
        fi (float): Incoming neutron flux correction factor.

    Returns:
        sigma_VWC (float): Uncertainty in terms of volumetric water content.

    References:
        Jakobi J, Huisman JA, Schrön M, Fiedler J, Brogi C, Vereecken H and Bogena HR (2020) Error Estimation for Soil Moisture Measurements With
        Cosmic Ray Neutron Sensing and Implications for Rover Surveys. Front. Water 2:10. doi: 10.3389/frwa.2020.00010
    """

    Ncorr = raw_counts * fw / (fp * fi)
    sigma_N = uncertainty_counts(raw_counts, metric="std", fp=fp, fw=fw, fi=fi)
    sigma_GWC = sigma_N * ((a0*N0) / ((Ncorr - a1*N0)**4)) * np.sqrt((Ncorr - a1 * N0)**4 + 8 * sigma_N**2 * (Ncorr - a1 * N0)**2 + 15 * sigma_N**4)
    sigma_VWC = sigma_GWC * bulk_density

    return sigma_VWC










